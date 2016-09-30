#!/usr/bin/python 
import utils
import numpy as np

if __name__ == "__main__":
    width = 50.0
    length = 5.0
    nav_status = np.array([10.0, 5.0, 1.57])
    print(utils.evaluate_sss_path(nav_status, width, length))
