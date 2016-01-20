#!/usr/bin/python

import pymoos
import time
import gmphd

comms = pymoos.comms()


def c():
    return comms.register('UHZ_HAZARD_REPORT', 0)


def m():
    map(lambda msg: msg.trace(), comms.fetch())
    return True


def main():

    comms.set_on_connect_callback(c)
    comms.set_on_mail_callback(m)
    comms.run('localhost', 9001, 'gmphd-moos')

    while True:
        time.sleep(1.0)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt as e:
        print("\n\tGMPHD-MOOS finished")
        pass
