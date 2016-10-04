#!/usr/bin/python 
import itertools
import utils
import numpy as np
import gmphd

if __name__ == "__main__":
    width = 50.0
    length = 5.0
    nav_status = np.array([10.0, 5.0, 1.57])
    print(utils.evaluate_sss_path(nav_status, width, length))

    # Test of the selection functions for features inside the fov of the sss 
    gmphd_components = [gmphd.GmphdComponent(1, [10,5],[1,0,0,1]),gmphd.GmphdComponent(1, [100,1],[1,0,0,1])]
    means = np.squeeze(np.asarray([comp._mean for comp in gmphd_components]))
    print(means,means.shape)
    sss_path = utils.evaluate_sss_path(nav_status, width, length)  
    gmm_mask = utils.inside_polygon(means, sss_path)
    gmm_masked = list(itertools.compress(gmphd_components, gmm_mask))
    print gmm_masked 
