import numpy as np
from scipy.stats import multivariate_normal
from copy import deepcopy
from operator import attrgetter
from itertools import *
# !/usr/bin/env python
# GM-PHD implementation  in Python by Dan Stowell updated by Tommaso Fabbri
#
# Based on the description in Vo and Ma (2006).
# (c) 2012 Dan Stowell and Queen Mary University of London.
# (c) 2016 Tommaso Fabbri and University of Pisa
# All rights reserved.
#
# NOTE: SPAWNING IS NOT IMPLEMENTED.

"""

This file is part of gmphd, GM-PHD filter in python by Dan Stowell.

    gmphd is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    gmphd is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with gmphd.  If not, see <http://www.gnu.org/licenses/>.
"""

# given an iterable of pairs return the key corresponding to the greatest value
def argmax(pairs):
    return max(pairs, key=itemgetter(1))[0]


# given an iterable of values return the index of the greatest value
def argmax_index(values):
    return argmax(enumerate(values))


# given an iterable of keys and a function f, return the key with largest f(key)
def argmax_f(keys, f):
    return max(keys, key=f)


class GmphdComponent:
    """

    GM-PHD Gaussian component.

    The Gaussian component is defined by:
        weight
        mean
        covariance

    """

    def __init__(self, weight, mean, cov):
        """
        """
        self.weight = np.float64(weight)

        self.mean = np.array(mean, dtype=np.float64, ndmin=2)
        self.cov = np.array(cov, dtype=np.float64, ndmin=2)

        self.mean = np.reshape(self.mean, (self.mean.size, 1))
        self.cov = np.reshape(self.cov, (self.mean.size, self.mean.size))


class GMPHD:
    birth_w = 0.001
    # birth_cov = np.array([[],[]]dtype=np.float64)
    """
    Represents a set of modelling parameters and the latest frame's
    GMM estimate, for a GM-PHD model without spawning.

    Typical usage would be, for each frame of input data, to run:
    g.update(obs)
    g.prune()
    estimate = g.extractstates()

    'gmm' is an array of GmphdComponent items which makes up
    the latest GMM, and updated by the update() call.
    It is initialised as empty.

    Test code example (1D data, with new trails expected at around 100):
    from gmphd import *
    g = Gmphd([GmphdComponent(1, [100], [[10]])], 0.9, 0.9, [[1]], [[1]], [[1]], [[1]], 0.000002)
    g.update([[30], [67.5]])
    g.gmmplot1d()
    g.prune()
    g.gmmplot1d()

    g.gmm

    [(float(comp.loc), comp.weight) for comp in g.gmm]
    """

    def __init__(self, birthgmm, survival, detection, f, q, h, r, clutter):
        """
            'gm' list of GmphdComponent

            'birthgmm' List of GmphdComponent items which makes up the GMM of birth probabilities.
            'survival' Survival probability.
            'detection' Detection probability.
            'f' State transition matrix F.
            'q' Process noise covariance Q.
            'h' Observation matrix H.
            'r' Observation noise covariance R.
            'clutter' Clutter intensity.
        """
        self.gm = []
        self.birthgmm = birthgmm

        self.survival = np.float64(survival)  # p_{s,k}(x) in paper
        self.detection = np.float64(detection)  # p_{d,k}(x) in paper

        self.f = np.array(f, dtype=np.float64)  # state transition matrix      (F_k-1 in paper)
        self.q = np.array(q, dtype=np.float64)  # process noise covariance     (Q_k-1 in paper)
        self.h = np.array(h, dtype=np.float64)  # observation matrix           (H_k in paper)
        self.r = np.array(r, dtype=np.float64)  # observation noise covariance (R_k in paper)
        self.clutter = np.float64(clutter)  # clutter intensity (KAU in paper)

    def create_birth(self, measures):
        born = [GmphdComponent(GMPHD.birth_w, m, self.r) for m in measures]
        return born

    def predict_birth(self, born_components):
        # Prediction for birth targets
        born = [GmphdComponent(comp.weight,
                               np.dot(self.f, comp.mean),
                               self.q + np.dot(np.dot(self.f, comp.cov), self.f.T)
                               ) for comp in born_components]
        return born

    def predict_existing(self):
        # Prediction for existing targets
        predicted = [GmphdComponent(self.survival * comp.weight,
                                    np.dot(self.f, comp.mean),
                                    self.q + np.dot(np.dot(self.f, comp.cov), self.f.T)
                                    ) for comp in self.gm]
        return predicted

    def update(self, measures, predicted):
        # Construction of PHD update components
        eta = [np.dot(self.h, comp.mean) for comp in predicted]
        s = [self.r + np.dot(np.dot(self.h, comp.cov), self.h.T) for comp in predicted]

        k = []
        for index, comp in enumerate(predicted):
            k.append(np.dot(np.dot(comp.cov, self.h.T), np.linalg.inv(s[index])))

        pkk = []
        for index, comp in enumerate(predicted):
            pkk.append(np.dot(np.eye(np.size(k[index])) - np.dot(k[index], self.h), comp.cov))

        # Update using the measures

        # The 'predicted' components are kept, with a decay
        pr_gm = [GmphdComponent(comp.weight * (1.0 - self.detection),
                                comp.mean, comp.cov) for comp in predicted]

        for z in measures:
            temp_gm = []
            for j, comp in enumerate(predicted):
                temp_gm.append(GmphdComponent(
                        self.detection * comp.weight * multivariate_normal(z, eta[j], s[j]),
                        comp.mean + np.dot(k[j], z - eta[j]),
                        comp.cov))

            # The Kappa thing (clutter and reweight)
            weight_sum = np.sum(comp.weight for comp in temp_gm)
            weight_factor = 1.0 / (self.clutter + weight_sum)
            for comp in temp_gm:
                comp.weight *= weight_factor
            pr_gm.extend(temp_gm)
        self.gm = pr_gm

    def run_iteration(self, measures, born_components):
        # Prediction for birthed targets
        pr_born = self.predict_birth(born_components)
        # Prediction for existing targets
        predicted = self.predict_existing()
        predicted.extend(pr_born)
        # Update
        self.update(measures, predicted)

    def prune(self, truncation_thresh=1e-6, merge_thresh=0.01, max_components=100):
        # Truncation step
        # weightsums = [simplesum(comp.weight for comp in self.gmm)]  # diagnostic
        truncated_gm = filter(lambda comp: comp.weight > truncation_thresh, self.gm)

        # weightsums.append(simplesum(comp.weight for comp in sourcegmm))
        # origlen = len(self.gmm)
        # trunclen = len(sourcegmm)

        # Iterate to build the new GMM
        newgmm = []
        while len(truncated_gm) > 0:
            # Find the component with the maximum weight
            j = argmax(comp.weight for comp in sourcegmm)
            comp_j = truncated_gm[j]
            truncated_gm.pop(j)
            # truncated_gm = truncated_gm[:j] + truncated_gm[j + 1:]

            # Computation of the Mahalanobis distances
            m_distances = []
            for i, comp in enumerate(truncated_gm):
                temp = np.dot((comp.mean - comp_j.mean).T, np.inv(comp.cov))
                m_distances.append(np.float64(np.dot(temp, (comp.mean - comp_j.mean))))
                if m_distances
            L = array([dist <= merge_thresh for dist in m_distances])
            subsumed = [weightiest]
            if any(dosubsume):
                # print "Subsuming the following locations into weightest with loc %s and weight %g (cov %s):" \
                #	% (','.join([str(x) for x in weightiest.loc.flat]), weightiest.weight, ','.join([str(x) for x in weightiest.cov.flat]))
                # print list([comp.loc[0][0] for comp in list(array(sourcegmm)[ dosubsume]) ])
                subsumed.extend(list(array(sourcegmm)[dosubsume]))
                sourcegmm = list(array(sourcegmm)[~dosubsume])
            # create unified new component from subsumed ones
            aggweight = simplesum(comp.weight for comp in subsumed)
            newcomp = GmphdComponent( \
                    aggweight,
                    sum(array([comp.weight * comp.loc for comp in subsumed]), 0) / aggweight,
                    sum(array([comp.weight * (comp.cov + (weightiest.loc - comp.loc) \
                                              * (weightiest.loc - comp.loc).T) for comp in subsumed]), 0) / aggweight
            )
            newgmm.append(newcomp)



        # Now ensure the number of components is within the limit, keeping the weightiest
        newgmm.sort(key=attrgetter('weight'))
        newgmm.reverse()
        self.gmm = newgmm[:maxcomponents]
        weightsums.append(simplesum(comp.weight for comp in newgmm))
        weightsums.append(simplesum(comp.weight for comp in self.gmm))
        print "prune(): %i -> %i -> %i -> %i" % (origlen, trunclen, len(newgmm), len(self.gmm))
        print "prune(): weightsums %g -> %g -> %g -> %g" % (weightsums[0], weightsums[1], weightsums[2], weightsums[3])
        # pruning should not alter the total weightsum (which relates to total num items) - so we renormalise
        weightnorm = weightsums[0] / weightsums[3]
        for comp in self.gmm:
            comp.weight *= weightnorm

    def extractstates(self, bias=1.0):
        """Extract the multiple-target states from the GMM.
		  Returns a list of target states; doesn't alter model state.
		  Based on Table 3 from Vo and Ma paper.
		  I added the 'bias' factor, by analogy with the other method below."""
        items = []
        print "weights:"
        print [round(comp.weight, 7) for comp in self.gmm]
        for comp in self.gmm:
            val = comp.weight * float(bias)
            if val > 0.5:
                for _ in range(int(round(val))):
                    items.append(deepcopy(comp.loc))
        for x in items: print x.T
        return items

    def extractstatesusingintegral(self, bias=1.0):
        """Extract states based on the expected number of states from the integral of the intensity.
		This is NOT in the GMPHD paper; added by Dan.
		"bias" is a multiplier for the est number of items.
		"""
        numtoadd = int(round(float(bias) * simplesum(comp.weight for comp in self.gmm)))
        print "bias is %g, numtoadd is %i" % (bias, numtoadd)
        items = []
        # A temporary list of peaks which will gradually be decimated as we steal from its highest peaks
        peaks = [{'loc': comp.loc, 'weight': comp.weight} for comp in self.gmm]
        while numtoadd > 0:
            windex = 0
            wsize = 0
            for which, peak in enumerate(peaks):
                if peak['weight'] > wsize:
                    windex = which
                    wsize = peak['weight']
            # add the winner
            items.append(deepcopy(peaks[windex]['loc']))
            peaks[windex]['weight'] -= 1.0
            numtoadd -= 1
        for x in items: print x.T
        return items

    ########################################################################################
    # def gmmeval(self, points, onlydims=None):
    #     """Evaluates the GMM at a supplied list of points (full dimensionality).
		# 'onlydims' if not nil, marginalises out (well, ignores) the nonlisted dims. All dims must still be listed in the points, so put zeroes in."""
    #     return [ \
    #         simplesum(comp.weight * comp.dmvnorm(p) for comp in self.gmm) \
    #         for p in points]
    #
    # def gmmeval1d(self, points, whichdim=0):
    #     "Evaluates the GMM at a supplied list of points (1D only)"
    #     return [ \
    #         simplesum(comp.weight * dmvnorm([comp.loc[whichdim]], [[comp.cov[whichdim][whichdim]]], p) for comp in
    #                   self.gmm) \
    #         for p in points]
    #
    # def gmmevalgrid1d(self, span=None, gridsize=200, whichdim=0):
    #     "Evaluates the GMM on a uniformly-space grid of points (1D only)"
    #     if span == None:
    #         locs = array([comp.loc[whichdim] for comp in self.gmm])
    #         span = (min(locs), max(locs))
    #     grid = (arange(gridsize, dtype=float) / (gridsize - 1)) * (span[1] - span[0]) + span[0]
    #     return self.gmmeval1d(grid, whichdim)
    #
    # def gmmevalalongline(self, span=None, gridsize=200, onlydims=None):
    #     """Evaluates the GMM on a uniformly-spaced line of points (i.e. a 1D line, though can be angled).
		# 'span' must be a list of (min, max) for each dimension, over which the line will iterate.
		# 'onlydims' if not nil, marginalises out (well, ignores) the nonlisted dims. All dims must still be listed in the spans, so put zeroes in."""
    #     if span == None:
    #         locs = array([comp.loc for comp in
    #                       self.gmm]).T  # note transpose - locs not a list of locations but a list of dimensions
    #         span = array(
    #                 [map(min, locs), map(max, locs)]).T  # note transpose - span is an array of (min, max) for each dim
    #     else:
    #         span = array(span)
    #     steps = (arange(gridsize, dtype=float) / (gridsize - 1))
    #     grid = array(map(lambda aspan: steps * (aspan[1] - aspan[0]) + aspan[0],
    #                      span)).T  # transpose back to list of state-space points
    #     return self.gmmeval(grid, onlydims)
    #
    # def gmmplot1d(self, gridsize=200, span=None, obsnmatrix=None):
    #     "Plots the GMM. Only works for 1D model."
    #     import matplotlib.pyplot as plt
    #     vals = self.gmmevalgrid1d(span, gridsize, obsnmatrix)
    #     fig = plt.figure()
    #     plt.plot(grid, vals, '-')
    #     fig.show()
    #     return fig
