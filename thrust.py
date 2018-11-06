import os
import sys
import numpy as np
import scipy.interpolate as ipl
from scipy import integrate

import matplotlib.pyplot as plt

#########################################################################
def thrust_load(motor_name):
    '''
    Load thrust data and interpolate thrust distribution

    input:
        motor_name : Load motor name
                     Defined in input file

    output:
        comb_time : combustion time [sec]
        linear_thrust : interpolated thrust distribution
                        (input:time, output:thrust force)
    '''

    # Load thrust data
    raw_data = np.loadtxt("input/thrust_data/"+motor_name, delimiter="\t", comments=";", skiprows=1)

    raw_time = raw_data[:,0]
    raw_thrust = raw_data[:,1]
    comb_time = raw_time[-1]

    # Interpolate distribution
    linear_thrust = ipl.interp1d(raw_time, raw_thrust, kind='linear',
                          bounds_error=False, fill_value=0.0)

    return comb_time, linear_thrust


#########################################################################
def show_hist(motor_name):
    '''
    Show thrust curve on graph
    '''

    comb_time, linear_thrust = thrust_load(motor_name)
    count = 1000
    time = np.linspace(0.0, comb_time, num=count)
    thrust = np.zeros(count)

    thrust_impulse = integrate.quad(linear_thrust, 0, comb_time)
    
    print("")
    #print("Max thrust    : {0}").format()
    print("Total impulse : {0} Ns".format(thrust_impulse[0]))

    for i in range(count):
        thrust[i] = linear_thrust(time[i])

    plt.figure()

    plt.plot(time,thrust)
    plt.xlabel("time[sec]")
    plt.ylabel("thrust[N]")
    plt.title(motor_name)
    plt.grid(color='grey')

    plt.show()


#########################################################################
if __name__ == "__main__":

    print("==> Select input mode")
    print("0: Use default setting (J250),  1:Use user setting")
    input_mode = input(">> ")
    input_mode = int(input_mode)

    if input_mode == 1:
        for num in range(5):
            print("")
            print("==> Input motor name from below list")
            files = os.listdir("./input/thrust_data/")
            print(files)
            print("")
            motor_name = input(">> ")

            if motor_name in files:
                print("")
                print("filename match!")
                break

            else:
                print("Input wrong file name. Please input again.")

    else:
        motor_name = 'Hypertek_J250.txt'

    show_hist(motor_name)

    