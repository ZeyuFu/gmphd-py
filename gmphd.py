import numpy as np
from operator import itemgetter, attrgetter
from scipy.stats import multivariate_normal
import math

# !/usr/bin/env python
# GM-PHD implementation  in Python by Dan Stowell modified by Tommaso Fabbri
#
# Based on the description in Vo and Ma (2006).
# (c) 2012 Dan Stowell and Queen Mary University of London.
# (c) 2016 Tommaso Fabbri and University of Pisa - Automation & Robotics Laboratory

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

np.set_printoptions(precision=3)


class GmphdComponent:
    """
    GM-PHD Gaussian component.

    The Gaussian component is defined by:
        _weight
        _mean
        _covariance
    """

    def __init__(self, weight, mean, cov):

        self._weight = np.float64(weight)
        self._mean = np.array(mean, dtype=np.float64)
        self._cov = np.array(cov, dtype=np.float64)

        print self._mean.size

        # self.mean.resize((self.mean.size, 1))
        # self.cov.resize((self.mean.size, self.mean.size))

    def __repr__(self):
        str_ = '\nW:\n{0}\nM:\n{1}\nC:\n{2}\n'.format(self._weight, self._mean, self._cov)
        return str_


class GMPHD:
    birth_w = 0.001

    def __str__(self):

        for i in self.gm:
            return i.__str__()

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

    def predict_birth(self, born_components):
        # Prediction for birth targets
        born = [GmphdComponent(comp._weight,
                               np.dot(self.f, comp._mean),
                               self.q + np.dot(np.dot(self.f, comp._cov), self.f.T)
                               ) for comp in born_components]
        return born

    def predict_existing(self):
        # Prediction for existing targets
        predicted = [GmphdComponent(self.survival * comp._weight,
                                    np.dot(self.f, comp._mean),
                                    self.q + np.dot(np.dot(self.f, comp._cov), self.f.T)
                                    ) for comp in self.gm]
        return predicted

    def update(self, measures, predicted):
        # Construction of PHD update components
        repr(predicted)

        eta = [np.dot(self.h, comp._mean) for comp in predicted]
        print 'Eta',  eta
        s = [self.r + np.dot(np.dot(self.h, comp._cov), self.h.T) for comp in predicted]

        k = []
        for index, comp in enumerate(predicted):
            k.append(np.dot(np.dot(comp._cov, self.h.T), np.linalg.inv(s[index])))

        pkk = []
        for index, comp in enumerate(predicted):
            pkk.append(np.dot(np.eye(np.shape(k[index])[0]) - np.dot(k[index], self.h), comp._cov))

        # Update using the measures

        # The 'predicted' components are kept, with a decay
        pr_gm = [GmphdComponent(comp._weight * (1.0 - self.detection),
                                comp._mean, comp._cov) for comp in predicted]

        for i in np.ndindex(measures.shape[1]):
            z = measures[:, i]
            temp_gm = []
            for j, comp in enumerate(predicted):
                print "Z", z.squeeze()
                print 'ETA', eta[j].squeeze()
                print 'S', s[j]
                mvn = multivariate_normal(eta[j].squeeze(), s[j])
                mvn_result = mvn.pdf(z.squeeze())
                temp_gm.append(GmphdComponent(
                        self.detection * comp._weight * mvn_result,
                        comp._mean + np.dot(k[j], z - eta[j]),
                        comp._cov))

            # The Kappa thing (clutter and reweight)
            weight_sum = np.sum(comp._weight for comp in temp_gm)
            if weight_sum != 0:
                weight_factor = 1.0 / (self.clutter + weight_sum)
                for comp in temp_gm:
                    comp._weight *= weight_factor
                pr_gm.extend(temp_gm)
        self.gm = pr_gm

    def run_iteration(self, measures, born_components):
        # Prediction for birthed targets
        print('Measures: ')
        print(measures)
        pr_born = self.predict_birth(born_components)
        # Prediction for existing targets
        predicted = self.predict_existing()
        predicted.extend(pr_born)
        print('Predicted components:'.format(predicted))
        # Update
        self.update(measures, predicted)
        print('Updated components:'.format(self.gm))
        # Prune
        self.prune()
        print('Pruning: '.format(self.gm))

    def prune(self, truncation_thresh=1e-6, merge_thresh=0.01, max_components=100):
        temp_sum_0 = np.sum([i._weight for i in self.gm])

        # Truncation step
        I = filter(lambda comp: comp._weight > truncation_thresh, self.gm)
        l = 0  # count the number of features/components
        pruned_gm = []

        # Merge step
        while len(I) > 0:
            l += 1
            j = np.argmax(i._weight for i in I)
            L = []
            indexes = []
            for index, i in enumerate(I):
                temp = np.dot((i._mean - I[j]._mean).T, np.linalg.inv(i._cov))
                mah_dist = np.float64(np.dot(temp, (i._mean - I[j]._mean)))
                if mah_dist <= merge_thresh:
                    L.append(i)
                    indexes.append(index)
            temp_weight = np.sum([i._weight for i in L])
            temp_mean = (1.0 / temp_weight) * np.sum([i._weight * i._mean for i in L], axis=0)
            temp_cov = np.zeros((temp_mean.size, temp_mean.size))

            for i in L:
                print 'TM', temp_mean
                print i._mean
                temp_cov += (i._cov + np.dot((temp_mean - i._mean).T, (temp_mean - i._mean)))
            pruned_gm.append(GmphdComponent(temp_weight, temp_mean, temp_cov))
            I = [i for j, i in enumerate(I) if j not in indexes]
        pruned_gm.sort(key=attrgetter('_weight'))
        pruned_gm.reverse()
        pruned_gm = pruned_gm[:max_components]
        temp_sum_1 = np.sum(i._weight for i in pruned_gm)
        for i in pruned_gm:
            i._weight *= temp_sum_0 / temp_sum_1

        self.gm = pruned_gm


def create_birth(measures):
    sigma_r = 2.0/3
    R = [[math.pow(2*sigma_r, 2), 0], [0, math.pow(2*sigma_r, 2)]]

    it = np.nditer(measures.shape[1])
    born = []
    for i in np.ndindex(measures.shape[1]):
        born.append(GmphdComponent(GMPHD.birth_w, measures[:,i], R))
    return born
