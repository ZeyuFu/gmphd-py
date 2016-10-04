#!/usr/bin/python

import pymoos
import json
import time
import math
import re
import numpy as np
import gmphd
import utils
from threading import Semaphore

measure = None
comms = None
f_gmphd = None
born_components = None
auv_nav_status = np.array([0,0,0], dtype=np.float64) # [x, y, yaw]

m_regex = re.compile('x=(?P<x>[-+]?(\d*[.])?\d+e?[-+]?\d*),y=(?P<y>[-+]?(\d*[.])?\d+e?[-+]?\d*),label=(?P<label>\d*),type=(?P<type>benign|hazard)')
msgs_semaphore = Semaphore()

def on_connect():
    global comms
    comms.register('NAV_X'            , 0)
    comms.register('NAV_Y'            , 0)
    comms.register('NAV_YAW'          , 0)
    comms.register('UHZ_HAZARD_REPORT', 0)
    return True

def on_mail():
    global auv_nav_status, measure
    msgs = comms.fetch()
    msgs_semaphore.acquire()
    for i in reversed(range(len(msgs))):
        if not msgs[i].string():
            if msgs[i].name() == 'UHZ_HAZARD_REPORT':
                g = m_regex.match(msgs[i].string())
                print('X={0},Y={1},L={2},T={3}'.format(g.group('x'), g.group('y'),
                                                       g.group('label'), g.group('type')))
                if g.group('type') == 'hazard':
                    measure = np.array([np.float64(g.group('x')), np.float64(g.group('y'))])
                    measure.shape = (measure.size, 1)
                    # on_measures(measure)
            if msgs[i].name() == 'NAV_X':
                auv_nav_status[0] = msgs[i].double()
            if msgs[i].name() == 'NAV_Y':
                auv_nav_status[1] = msgs[i].double()
            if msgs[i].name() == 'NAV_YAW':
                auv_nav_status[2] = msgs[i].double()
    msgs_semaphore.release()
    return True


def on_measures(measure):
    # Measures is an array of [x, y] of possible
    # TODO Implement the call to the filter
    global f_gmphd, born_components

    f_gmphd.run_iteration(measure, born_components)
    # TODO Feature all'interno del FOV hanno una evoluzione differente
    # f_gmphd = f_gmphd.run_iteration(measure, born_components)

    born_components = gmphd.create_birth(measure)
    print('Born components: {0}'.format(born_components))


def main(_sigma_q, _sigma_r, _p_d, _p_s):
    global comms, f_gmphd, born_components
    comms = pymoos.comms()

    F = [[1, 0], [0, 1]]
    H = [[1, 0], [0, 1]]

    sigma_q = _sigma_q
    sigma_r = _sigma_r

    Q = [[math.pow(sigma_q, 2), 0], [0, math.pow(sigma_q, 2)]]
    R = [[math.pow(2*sigma_r, 2), 0], [0, math.pow(2*sigma_r, 2)]]

    p_d = _p_d
    p_s = _p_s

    clutter_intensity = 0

    born_components = []

    f_gmphd = gmphd.GMPHD([], p_s, p_d, F, Q, H, R, clutter_intensity)

    comms.set_on_connect_callback(on_connect)
    comms.set_on_mail_callback(on_mail)
    comms.run('localhost', 9001, 'gmphd-moos')
    i = 0
    while True:
        i = i + 1
        # print(i)
        time.sleep(.01)

if __name__ == "__main__":
    with open('filter-conf.json') as conf_file:
        data = json.load(conf_file)
    sigma_q = data['Filter']['sigma_q']
    sigma_r = data['Filter']['sigma_r']
    p_d = data['Filter']['prob_d']
    p_s = data['Filter']['prob_s']

    try:
        main(sigma_q, sigma_r, p_d, p_s)
    except KeyboardInterrupt as e:
        print("\nGMPHD-MOOS terminated")
        pass
