#!/usr/bin/env python

"""
.. module:: runtime
    :synopsis: Tools to gather info about runtime enviroment,
               currently contains only the number of CPUs

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def nCPUs():
    """ obtain the number of CPU cores on the machine, for several
        platforms and python versions. """
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except Exception,e:
        pass
    try:
        import psutil
        return psutil.NUM_CPUS
    except Exception,e:
        pass
    try:
        import os
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))
        if res>0: return res
    except Exception,e:
        pass
    return None


if __name__ == "__main__":
    print "This machine has %d CPUs" % nCPUs()
