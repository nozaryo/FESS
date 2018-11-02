"""
F.T.E. Rocket Simulating System
Edited by Ryoya Nozaki

reference
https://github.com/Shun0228/RocketSimulation
https://qiita.com/taashi/items/400871fb13df476f42d2

"""
import os
import sys
import math
from argparse import ArgumentParser
import numpy as np
import scipy as sp
import numpy.linalg as nplin
from scipy.integrate import odeint

import calculation as cal
import quaternion as qt
import environment as env
import summary as summ


def get_option():
    args = ArgumentParser(prog='F.T.E. Rocket Simulating System')
    args.add_argument('-J','--json',default="",type=str,
                           help='json file name ')
    args.add_argument('-M','--mood',
                            type=str, choices=['d','s'], nargs=1, metavar='mood',
                            help="Select mood d:detail s:scatter")
    args.add_argument('-P','--parachute',
                            type=int, choices=[0,1], nargs=1, metavar='open_mood',
                            help="Select parachute mood 0: Trajectory,  1: Open one stage parachute")
    return args.parse_args()

if __name__ == '__main__':
    args = get_option()
    print(args)
