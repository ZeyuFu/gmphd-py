#!/usr/bin/python

import pymoos
import time
import math
import re
import numpy as np
import gmphd

comms = None
f_gmphd = None
born_components = None


def on_connect():
    global comms
    return comms.register('UHZ_HAZARD_REPORT', 0)


def on_mail():
    m_regex = re.compile('x=(?P<x>[-+]?(\d*[.])?\d+e?[-+]?\d*),'
                         'y=(?P<y>[-+]?(\d*[.])?\d+e?[-+]?\d*),label=(?P<label>\d*),'
                         'type=(?P<type>benign|hazard)')

    msgs = comms.fetch()
    for i in reversed(range(len(msgs))):
        if msgs[i].name() == 'UHZ_HAZARD_REPORT':
            g = m_regex.match(msgs[i].string())
            print('X={0},Y={1},L={2},T={3}'.format(g.group('x'), g.group('y'),
                                                   g.group('label'), g.group('type')))
            if g.group('type') == 'hazard':
                measure = np.array([np.float64(g.group('x')), np.float64(g.group('y'))])
                measure.shape = (measure.size, 1)
                print(measure)
                on_measures(measure)
    return True


def on_measures(measure):
    # Measures is an array of [x, y] of possible
    # TODO Implement the call to the filter
    global f_gmphd, born_components
    print(measure)
    f_gmphd.run_iteration(measure, born_components)
    # f_gmphd = f_gmphd.run_iteration(measure, born_components)

    born_components = gmphd.create_birth(measure)
    print('Born components: '.format(born_components))


def main():
    global comms, f_gmphd, born_components
    comms = pymoos.comms()

    # TODO Implement a configuration file *.csv for speed up the simulations and settings

    F = [[1, 0], [0, 1]]
    sigma_q = 1e-3
    Q = [[math.pow(sigma_q, 2), 0], [0, math.pow(sigma_q, 2)]]

    H = [[1, 0], [0, 1]]
    sigma_r = 2.0/3
    R = [[math.pow(2*sigma_r, 2), 0], [0, math.pow(2*sigma_r, 2)]]

    p_d = 0.95
    p_s = 1.00

    clutter_intensity = 0

    born_components = []

    f_gmphd = gmphd.GMPHD([], p_s, p_d, F, Q, H, R, clutter_intensity)

    comms.set_on_connect_callback(on_connect)
    comms.set_on_mail_callback(on_mail)
    comms.run('localhost', 9001, 'gmphd-moos')

    while True:
        time.sleep(1.0)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt as e:
        print("\n\tGMPHD-MOOS terminated")
        pass
