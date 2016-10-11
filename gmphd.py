#!/usr/bin/env python

import numpy as np
from operator import itemgetter, attrgetter
from scipy.stats import multivariate_normal
import math
from matplotlib import path
import utils
import itertools
import copy
# GM-PHD implementation in Python by Dan Stowell modified by Tommaso Fabbri
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
        self._mean.resize((self._mean.size, 1))
        self._cov.resize((self._mean.size, self._mean.size))

    def __repr__(self):
        str_ = 'W: {0} M: {1} C: {2}'.format(self._weight, self._mean.tolist(), self._cov.tolist())
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

        self.clutter = np.float64(clutter)  # clutter intensity
        self.swath_width = np.float64(50)   # uFldSensor Swath Width
        self.swath_length = np.float64(4)   # uFldSensor Swath Length

    def update_ufldsensor_data(self, width, length):
        self.swath_width = width
        self.swath_length = length

    def predict_birth(self, born_components):
        # Prediction for birth targets
        born = [GmphdComponent(comp._weight,
                               np.dot(self.f, comp._mean),
                               self.q + np.dot(np.dot(self.f, comp._cov), self.f.T)
                               ) for comp in born_components]
        return born

    def predict_existing(self, nav_status, sss_path):
        # Prediction for existing targets
        if len(self.gm) == 0:
            return []
        means = np.asarray(np.array([ comp._mean for comp in self.gm]))
        means = np.squeeze(means)
        means = np.squeeze(means)
        means = np.reshape(means, (means.size / 2, 2))
        if len(means) > 0:
            gmm_mask = utils.inside_polygon(means, sss_path)
            gmm_fov = list(itertools.compress(self.gm, gmm_mask))
            predicted = []
            for idx, comp in enumerate(self.gm):
                if gmm_mask[idx]:
                    predicted.append(GmphdComponent(self.survival * comp._weight,
                                               np.dot(self.f, comp._mean), self.q + np.dot(np.dot(self.f, comp._cov), self.f.T))) 
                else:
                    predicted.append(self.gm[idx])
        return predicted

    def auv_update(self, measures, nav_status, sss_path):
        '''
            Construction of the update components
            s   = R + H * _Cov * H.T 
            eta = H * _Mean
            K = _Cov * H.T * ( H * _Cov * H.T + R) ^ -1 
            Pkk = ( I - K * H ) * _Cov
        '''
        if len(self.gm) == 0:
           return 

        eta = [np.dot(self.h, comp._mean) for comp in self.gm]
        s = [self.r + np.dot(np.dot(self.h, comp._cov), self.h.T) for comp in self.gm]

        k = [] 
        for index, comp in enumerate(self.gm):
            k.append(np.dot(np.dot(comp._cov, self.h.T), np.linalg.inv(s[index])))

        pkk = []
        for index, comp in enumerate(self.gm):
            pkk.append(np.dot(np.eye(np.shape(k[index])[0]) - np.dot(k[index], self.h), comp._cov))

        '''
            The predicted components are kept with a decay (1 - Pd) iff they are inside the FOV, contained into pr_gm
        '''

        # Extraction of the coordinates of the components inside the RFS 
        means = np.asarray(np.array([ comp._mean for comp in self.gm]))
        means = np.squeeze(means)
        means = np.reshape(means, (means.size / 2, 2))
        # Temp copy of the RFS after the prediction step
        pr_gm = copy.deepcopy(self.gm) 
        # Mask of the components inside the FOV
        gmm_mask = utils.inside_polygon(means, sss_path)
        
        # print(gmm_mask.shape, pr_gm)  

        # gmm_fov = list(itertools.compress(self.gmm, gmm_mask))
        for idx, comp in enumerate(pr_gm):
            if gmm_mask[idx]:
               pr_gm[idx]._weight *= (1 - self.detection) 

        for i in np.ndindex(measures.shape[1]):
            z = measures[:, i]
            temp_gm = []
            for j, comp in enumerate(self.gm):
                # Computation of q_k
                mvn = multivariate_normal(eta[j].squeeze(), s[j])
                mvn_result = mvn.pdf(z.squeeze())

                temp_gm.append(GmphdComponent(
                        self.detection * comp._weight * mvn_result,
                        comp._mean + np.dot(k[j], z - eta[j]),
                        comp._cov))

            # The Kappa thing (clutter and reweight)
            weight_sum = np.sum(comp._weight for comp in temp_gm)
            
            if weight_sum >= 1e-9:
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

    def run_auv_iteration(self, measures, born_components, nav_status, swath_w, swath_l):
        # pr_born = self.predict_birth(born_components)
        sss_path = utils.evaluate_sss_path(nav_status, swath_w, swath_l)  
      
        predicted = self.predict_existing(nav_status, sss_path)
        predicted.extend(born_components)
        self.gm = predicted

        self.auv_update(measures, nav_status, sss_path)
        self.prune()


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
            # Loop for compute the Mahalanobis distance among the components survived the trunctation step 
            for index, i in enumerate(I):
                temp = np.dot((i._mean - I[j]._mean).T, np.linalg.inv(i._cov))
                mah_dist = np.float64(np.dot(temp, (i._mean - I[j]._mean)))
                if mah_dist <= merge_thresh:
                    L.append(i)
                    indexes.append(index)
            # L contains the components with Mahalanobis distance greater than merge_thresh 
            # between component j and the others         
            temp_weight = np.sum([i._weight for i in L])
            temp_mean = (1.0 / temp_weight) * np.sum([i._weight * i._mean for i in L], axis=0)
            temp_cov = np.zeros((temp_mean.size, temp_mean.size))
            for i in L:
                # print 'TM', temp_mean.tolist()
                # print i._mean.tolist()
                temp_cov += (i._cov + np.dot((temp_mean - i._mean).T, (temp_mean - i._mean)))
                
            pruned_gm.append(GmphdComponent(temp_weight, temp_mean, temp_cov))
            I = [i for j, i in enumerate(I) if j not in indexes]
        pruned_gm.sort(key=attrgetter('_weight'))
        pruned_gm.reverse()
        pruned_gm = pruned_gm[:max_components]
        # temp_sum_1 = np.sum(i._weight for i in pruned_gm)
        # for i in pruned_gm:
        #     i._weight *= temp_sum_0 / temp_sum_1
        self.gm = pruned_gm


def create_birth(measures):
    sigma_r = 2.0/3
    R = [[math.pow(2*sigma_r, 2), 0], [0, math.pow(2*sigma_r, 2)]]

    it = np.nditer(measures.shape[1])
    born = []
    for i in np.ndindex(measures.shape[1]):
        born.append(GmphdComponent(GMPHD.birth_w, measures[:,i], R))
    return born

