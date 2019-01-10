"""
F.T.E. Rocket Simulating System
Edited by Ryoya Nozaki

reference
https://github.com/Shun0228/RocketSimulation

"""
import os
import sys
import math
import json
import datetime
from argparse import ArgumentParser
import numpy as np
import scipy as sp
import numpy.linalg as nplin
from scipy.integrate import odeint

import calculation as cal
import quaternion as qt
import environment as env
import plot_process as plot

# Define computation setting ---------------------------
dt = 0.001        # time interval
end_time = 2000.0        # computation end time

# Define wind pattern setting
wind_vel_st = 1.0    # minimum wind velocity [m/s]
wind_vel_interval = 1.0    # interval of wind velocity [m/s]
vel_pat = 7   # velocity range [m/s]
dir_pat = 8   # wind direction derivation (deg = 360.0/dir_pat)

# Define launcher elevation setting
elev_st = 80.0    # Launcher elevation start angle [deg]
elev_interval = 1.0    # Launche elevation angle interval [deg]
elev_pat = 3    # Launcher elevation range

# Define vector ---------------------------------------
time_vec = np.arange(0, end_time, dt)
case_total = vel_pat * dir_pat
drop_point = np.zeros([2, dir_pat+1, vel_pat])
wind_case = np.array([wind_vel_st, wind_vel_interval, vel_pat, dir_pat])

# Generate rsim object ----------------------------------------
sim = cal.RocketSim()
post =  plot.PlotProcess()


def get_option():
    args = ArgumentParser(prog='F.T.E. Rocket Simulating System')

    args.add_argument('-J','--json',default="",type=str, metavar='file_name',
                           help='json file name ')

    args.add_argument('-M','--mode',
                            type=str, choices=['d','s','sm','ss'], nargs=1, metavar='mode',
                            help="Select mode d:detail s:scatter ss:scatter single sm:scatter multiple")

    args.add_argument('-P','--parachute',
                            type=int, choices=[0,1], nargs=1, metavar='open_mood',
                            help="Select parachute mode 0: Trajectory,  1: Open one stage parachute")

    args.add_argument('-W','--wind',
                            type=float, nargs=2, metavar=('vel','dir'),
                            help="input wind velocity and direction")

    return args.parse_args()

if __name__ == '__main__':
    args = get_option()

files = os.listdir("./input/")

#if '.json' in args.json[0]:
#    args.json[0] = args.json[0]+".json"

if not args.json in files:

    if args.json:
        print("")
        print("Input wrong file name. Please input again.")

    filenames = list()
    for file in files:
            if '.json' in file:
                filenames.append(file.replace('.json',''))

    print("")
    print("==> Select rocket from below lists, and input rocket name")
    print(filenames)

    for num in range(5):
        filename = input('>> ')
        if '.json' in filename:
            if filename in files:
                break
        else:
            filename = filename + '.json'
            if filename in files:
                break
    else:
        print("Input wrong rocket name. Please input again.")
else:
    filename = args.json

try:
    # Load parameter from user setting
    f = open('input/'+filename, 'r')
    stdin_json = json.load(f)
    f.close()

    stdin_info = stdin_json["info"]
    stdin_rocket = stdin_json["rocket"]
    stdin_motor = stdin_json["motor"]
    stdin_env = stdin_json["environment"]

except json.decoder.JSONDecodeError:
    print("===========================================================================")
    print("[ERROR]: Json Decode Error")
    print("Lack of delimiter ','")
    print("Check ./input/****.json")
    print("===========================================================================")

except KeyError:
    print("===========================================================================")
    print("[ERROR]: Wrong dictionary name")
    print("There is incorrect name.")
    print("Please check dictionary name at ./input/****.json.")
    print("'rocket','motor' and 'environment' is correct name")
    print("==========================================================================")
    sys.exit()

# Write input infomation
print("")
print("Team    : " + stdin_info.get("TEAM"))
print("Name    : " + stdin_info.get("NAME"))
print("Date    : " + str(datetime.date.today()))
print("Version : " + stdin_info.get("VERSION"))
print("")
print("Motor   : " + stdin_motor["motor_file"])
print("")
print("Place   : " + stdin_env.get("place"))
print("rail_len: " + str(stdin_env.get("rail_len")))
print("rail_azi: " + str(stdin_env.get("rail_azi")))
print("")
# Load rocket parameter
try:
    print("ref_len : " + str(stdin_rocket["ref_len"]))
    print("diam    : " + str(stdin_rocket["diam"]))
    print("CGlen_i : " + str(stdin_rocket["CGlen_i"]))
    print("CGlen_f : " + str(stdin_rocket["CGlen_f"]))
    print("mass_i  : " + str(stdin_rocket["mass_i"]))
    print("mass_f  : " + str(stdin_rocket["mass_f"]))
    print("Iyz_i   : " + str(stdin_rocket["Iyz_i"]))
    print("Iyz_f   : " + str(stdin_rocket["Iyz_f"]))
    print("CPlen   : " + str(stdin_rocket["CPlen"]))
    print("Cd      : " + str(stdin_rocket["Cd"]))
    print("Cna     : " + str(stdin_rocket["Cna"]))


except KeyError:
    print("===========================================================================")
    print("[ERROR]: lack of parameters")
    print("There is no data of rocket.")
    print("Please check parameter at ./input/****.json.")
    print("==========================================================================")
    sys.exit()


if "Cmq" not in stdin_rocket:
    Cmq = -(stdin_rocket["Cna"]/2)*( ( (stdin_rocket["CPlen"] - stdin_rocket["CGlen_f"])/ stdin_rocket["ref_len"] ) )**2
    print("Cmq     : " + str(Cmq))
else:
    print("Cmq     : " + str(stdin_rocket.get("Cmq", -2.0)))

# Load 1st parachute parameter
if "vel_1st" not in stdin_rocket:
    print("---------------------------------------------------------------------------")
    print("[WARNING]: lack of parameters")
    print("There is no data of 1st parachute.")
    print("If you want to calculate one stage parachute, please add below parameter.")
    print("vel_1st, op_type_1st, op_time_1st, op_height_1st")
    print("---------------------------------------------------------------------------")

print("vel_1st : " + str(stdin_rocket["vel_1st"]))
print("delay_time_1st: " + str(stdin_rocket.get("delay_time_1st", 0.0)))
print("vel_2nd : " + str(stdin_rocket.get("vel_2nd", 0.0)))    # 0 means not to open 2nd para

sim.set_param(filename)
if not args.mode:
    print("")
    print("==> Select Simulation Mode")
    print("0: Detail Result,  1: Scatter Result")
    sim_mode = input('>> ')
    sim_mode = int(sim_mode)
elif 'd' in args.mode[0]:
    sim_mode = 0
elif 's' in args.mode[0]:
    sim_mode = 1


if not args.parachute:
    print("")
    print("==> Select parachute mode")
    print("0: Trajectory,  1: Open parachute")
    op_flg = input('>> ')
    op_flg = int(op_flg)
else:
    op_flg = args.parachute[0]

if sim_mode == 1:
    if not args.mode or args.mode[0] == 's':
        print("")
        print("==> Select Elevation Mode")
        print("0: Single output,  1: Multiple output")
        elev_mode = input('>> ')
        elev_mode = int(elev_mode)
    elif 'ss' in args.mode[0]:
        elev_mode = 0
    elif 'sm' in args.mode[0]:
        elev_mode = 1

else:
    elev_mode = 0

if elev_mode == 0:
    if not "rail_elev" in stdin_env:
        print("")
        print("==> Select Launcher Elevation (Vertical:90deg)")
        rail_elev = input('>> ')
        rail_elev = float(rail_elev)
    else:
        rail_elev = float(stdin_env.get("rail_elev"))
    elev_pat = 1
elif elev_mode == 1:
    rail_elev = elev_st
    print("")
    print("Launcher elevation : {0} - {1} deg".format(rail_elev, rail_elev + elev_pat -1.0))

rail_len = float(stdin_env.get("rail_len"))

rail_azi = float(stdin_env.get("rail_azi"))

rail_cond = np.array([rail_len, rail_elev, rail_azi])

# Execute Rocket simulation ------------------------------------
if sim_mode == 0:
    # Detail mode

    # Set wind condition
    if not args.wind:
        print("")
        print("==> Input wind velocity [m/s]")
        wind_vel = input('>> ')
        wind_vel = float(wind_vel)

        print("")
        print("==> Input wind direction [deg] (East:0deg, North:90deg)")
        wind_dir = input('>> ')
        wind_dir = float(wind_dir)
    else:
        wind_vel = args.wind[0]
        wind_dir = args.wind[1]

    wind_cond = np.array([wind_vel, wind_dir])

    # Compute ode equation
    sim.set_environment(rail_cond, wind_cond)
    sim.set_orbit_condition(op_flg)
    result = sim.execute(time_vec)

    # Post process
    post.set_variety(result, wind_cond,rail_cond)
    post.plot_detail()

elif sim_mode == 1:
    # Scatter plot mode

    print("")
    print("Calculation start")

    for ie in range(elev_pat):

        rail_cond[1] = rail_elev + ie

        print("")
        print("--- Elevation = {0} deg --------------------".format(rail_cond[1]))

        for iv in range(vel_pat):

            print("--------------------------------------------")

            for id in range(dir_pat):

                wind_vel = wind_vel_st + iv * wind_vel_interval
                wind_dir = id * (360.0/dir_pat)
                wind_cond = np.array([wind_vel, wind_dir])
                print('Wind vel:{0[0]:>4} m/s,   Wind dir:{0[1]:>6} deg'.format(wind_cond))

                # Compute ode equation
                sim.set_environment(rail_cond, wind_cond)
                sim.set_orbit_condition(op_flg)
                result = sim.execute(time_vec)

                # Write place of landing point
                drop_point[:, id, iv] = result[-1, 4:6]

                # Store result
                post.set_variety(result, wind_cond,rail_cond)


        # Copy overlapped value for plotting
        drop_point[:, -1, :] = drop_point[:, 0, :]

        for i in range(int(vel_pat/wind_vel_interval)) :
            for j in range(dir_pat) :
                print([i,j,drop_point[0,j,i],drop_point[1,j,i]])

        # Post process
        post.set_map(stdin_env.get("place"),rail_cond)
        post.plot_scatter(filename,wind_case,rail_cond[1],ie,op_flg,elev_mode)

    print("")
    print("Calculation end!")
