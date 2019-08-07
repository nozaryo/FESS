import sys
import math
import json
import datetime
import numpy as np
import scipy as sp
import numpy.linalg as nplin
import quaternion as qt
import thrust as th
import environment
from scipy.integrate import odeint

#########################################################################
class RocketSim:
    """ 5DoF rocket simulation """

    #########################################################################
    def set_param(self, filename):
        """ Overwrite rocket parameter when this method is called """
        self.env = environment.setEnv()
        self.env.wind_setting()


        self.peak_time = 0
        try:
            # Load parameter from user setting
            f = open('input/'+filename, 'r')
            stdin_json = json.load(f)
            f.close()

            stdin_info = stdin_json["info"]
            stdin_rocket = stdin_json["rocket"]
            stdin_motor = stdin_json["motor"]
            stdin_env = stdin_json["environment"]
            self.windFileName = str(stdin_env.get("wind_file", 0));
            self.env.wind_file_set(str(stdin_env.get("wind_file", 0)))

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
            print("'rocket' and 'motor' is correct name")
            print("==========================================================================")
            sys.exit()

        # Load rocket parameter
        self.Ix_i = 0.02
        self.Ix_f = 0.01
        try:
            self.ref_len = stdin_rocket["ref_len"]
            self.diam = stdin_rocket["diam"]
            self.CGlen_i = stdin_rocket["CGlen_i"]
            self.CGlen_f = stdin_rocket["CGlen_f"]
            self.mass_i = stdin_rocket["mass_i"]
            self.mass_f = stdin_rocket["mass_f"]
            self.Iyz_i = stdin_rocket["Iyz_i"]
            self.Iyz_f = stdin_rocket["Iyz_f"]
            self.CPlen = stdin_rocket["CPlen"]
            self.Cd = stdin_rocket["Cd"]
            self.Cna = stdin_rocket["Cna"]


        except KeyError:
            print("===========================================================================")
            print("[ERROR]: lack of parameters")
            print("There is no data of rocket.")
            print("Please check parameter at ./input/****.json.")
            print("==========================================================================")
            sys.exit()

        if "Cmq" not in stdin_rocket:
            self.Cmq = -(stdin_rocket["Cna"]/2.0)*( ( (stdin_rocket["CPlen"] - stdin_rocket["CGlen_f"])/ stdin_rocket["ref_len"] ) )**2
        else:
            self.Cmq = stdin_rocket.get("Cmq", -2.0)


        # Load 1st parachute parameter
        if "vel_1st" not in stdin_rocket:
            print("---------------------------------------------------------------------------")
            print("[WARNING]: lack of parameters")
            print("There is no data of 1st parachute.")
            print("If you want to calculate one stage parachute, please add below parameter.")
            print("vel_1st, op_type_1st, op_time_1st, op_height_1st")
            print("---------------------------------------------------------------------------")

        self.vel_1st = stdin_rocket["vel_1st"]
        self.op_type_1st = stdin_rocket["op_type_1st"]
        self.op_time_1st = stdin_rocket.get("op_time_1st")
        self.delay_time_1st = stdin_rocket.get("delay_time_1st", 0.0)


        # Load 2nd parachute parameter
        if "vel_2nd" not in stdin_rocket:
            print("---------------------------------------------------------------------------")
            print("[WARNING]: lack of parameters")
            print("There is no data of 2nd parachute.")
            print("If you want to calculate two stage parachute, please add below parameter.")
            print("vel_2nd, op_type_2nd, op_time_2nd, op_height_2nd")
            print("---------------------------------------------------------------------------")

        self.vel_2nd = stdin_rocket.get("vel_2nd", 0.0)   # 0 means not to open 2nd para
        self.op_type_2nd = stdin_rocket.get("op_type_2nd")

        if self.op_type_2nd == 0:
            try:
                self.op_height_2nd = stdin_rocket["op_height_2nd"]
            except KeyError:
                print("===========================================================================")
                print("[ERROR]: Lack of parameter")
                print("There is no data of 2nd para open height")
                print("Add op_height_2nd in json file")
                print("===========================================================================")
                sys.exit()

        elif self.op_type_2nd == 1:
            try:
                self.op_time_2nd = stdin_rocket["op_time_2nd"]
            except KeyError:
                print("===========================================================================")
                print("[ERROR]: Lack of parameter")
                print("There is no data of 2nd para open time")
                print("Add op_time_2nd in json file")
                print("===========================================================================")
                sys.exit()

        # Load 3rd parachute parameter
        self.vel_3rd = stdin_rocket.get("vel_3rd", 0.0)    # 0 means not to open 3rd para
        self.op_type_3rd = stdin_rocket.get("op_type_3rd")

        if self.op_type_3rd == 0:
            try:
                self.op_height_3rd = stdin_rocket["op_height_3rd"]
            except KeyError:
                print("===========================================================================")
                print("[ERROR]: Lack of parameter")
                print("There is no data of 3rd para open height")
                print("Add op_height_3rd in json file")
                print("===========================================================================")
                sys.exit()

        elif self.op_type_3rd == 1:
            try:
                self.op_time_3rd = stdin_rocket["op_time_3rd"]
            except KeyError:
                print("===========================================================================")
                print("[ERROR]: Lack of parameter")
                print("There is no data of 3rd para open time")
                print("Add op_time_3rd in json file")
                print("===========================================================================")
                sys.exit()

        # Set 4th
        self.op_time_4th = stdin_rocket.get("op_time_4th")
        self.op_height_4th = stdin_rocket.get("op_height_4th")


        # Load motor parameter
        self.motor_name = stdin_motor["motor_file"]
        self.motor_2nd = stdin_motor.get("motor_2nd")
        self.iject_time = stdin_motor.get("iject_time_2nd")
        self.iject_height = stdin_motor.get("iject_height_2nd")

        # Load thrust data
        self.comb_time, self.linear_thrust = th.thrust_load(self.motor_name)

        if self.motor_2nd is not None:
            self.comb_time_2nd, self.linear_thrust_2nd = th.thrust_load(self.motor_2nd)

        # Compute other parameters from above ones
        self.ref_sq = math.pi * 0.25 * self.diam ** 2
        self.para_sq = 1.0
        self.Cd_para_1st = self.mass_f * 9.81 / (0.5*1.25*(self.vel_1st**2)*self.para_sq)

        if self.vel_2nd > 0.0:
            self.Cd_para_2nd = self.mass_f * 9.81 / (0.5*1.25*(self.vel_2nd**2)*self.para_sq)
        else:
            self.Cd_para_2nd = 0.0


        if self.vel_3rd > 0.0:
            self.Cd_para_3rd = self.mass_f * 9.81 / (0.5*1.25*(self.vel_3rd**2)*self.para_sq)
        else:
            self.Cd_para_3rd = 0.0

        self.max_height = 0.0




    #########################################################################
    def set_environment(self, launcher_condition, wind_condition):
        """ Overwrite environment parameter when this method is called """

        # Set launcher parameter
        self.rail_len = launcher_condition[0]
        self.rail_elev = launcher_condition[1]
        self.rail_azi = launcher_condition[2]

        # Set wind parameter
        self.wind_ref = wind_condition[0]
        self.wind_dir = wind_condition[1]
        self.wind = np.array([self.wind_ref, self.wind_dir])

    #########################################################################
    def set_orbit_condition(self, open_condition):
        """ Overwrite open condition parameter when this method is called """

        # Set launcher parameter
        self.op_flg = open_condition

    #########################################################################
    def rocket_dynamics(self, var, time):

        """
        Function of calculating rocket dynamics only for trajectory orbit

        Parameter:
            time : time [sec]

            var  : variety vector
                    var[0] : mass [kg]
                    var[1] : length from nose to center of mass [m]
                    var[2] : inertia moment of pitching & yawing [kg*m^2]
                    var[3] : inertia moment of rolling [kg*m^2]
                    var[4:7] : position [m] (ENU coordinate)
                    var[7:10] : velocity [m/s]
                    var[10:13] : angular velocity (roll,pitch,yaw)
                    var[13:17] : quaternion

            wind : wind vector
                    wind[0]   : wind velocity [m/s]
                    wind[1]   : wind direction [deg]

        Return:
            dx : time derivative of variety

        """
        # Set each function from var
        self.env.wind_file_set(self.windFileName)
        mass = var[0]
        dif_len = self.CPlen - var[1]
        Iyz = var[2]
        Ix = var[3]
        pos = var[4:7]
        vel = var[7:10]
        omega = var[10:13]
        quat = qt.quat_normalize(var[13:17])

        # Define initial vector
        dVar = np.empty(var.size)
        force_b = np.empty(3)
        mom_damp = np.empty(3)
        mom_b = np.empty(3)
        dcm = qt.quat2dcm(quat)
        Gravity = np.array([0.0, 0.0, -9.81])
        wind_vel = self.env.wind_method(pos[2], self.wind)
        temp, pres, dens = self.env.gas_param(pos[2])

        # Define velocity vector and norm
        vel_b = dcm @ (vel - wind_vel)
        vel_norm = nplin.norm(vel_b)

        # Define local angle ---------------------------------------------
        if vel_norm == 0:
            alpha = 0.0
            beta = 0.0
            gamma = 0.0

        else:
            alpha = np.arctan(vel_b[2]/vel_b[0])
            beta = np.arctan(vel_b[1]/vel_b[0])
            gamma = 0.0


        # Rocket parameter condition --------------------------------------
        if time <= self.comb_time:
            # Decrease mass by combustion
            dVar[0] = (self.mass_f-self.mass_i) / self.comb_time
            dVar[1] = (self.CGlen_f-self.CGlen_i) / self.comb_time
            dVar[2] = (self.Iyz_f-self.Iyz_i) / self.comb_time
            dVar[3] = (self.Ix_f-self.Ix_i) / self.comb_time

        else:
            # Finish combustion
            dVar[0] = 0.0
            dVar[1] = 0.0
            dVar[2] = 0.0
            dVar[3] = 0.0


        # Aerodynamics coefficient condition ------------------------------
        if nplin.norm(pos) <= self.rail_len and vel[2] >= 0.0:
            # On launcher rail
            Cf_rail = 0.05
            Cd = self.Cd + Cf_rail

            Cn_p = 0.0
            Cn_y = 0.0
            Cmq_p = 0.0
            Cmq_y = 0.0

        # Store trajectory peak time
        if pos[2] > self.max_height:
            self.max_height =  pos[2]
            self.peak_time = time

        # Coasting
        Cd = self.Cd
        Cn_p = self.Cna * alpha
        Cn_y = self.Cna * beta
        Cmq_p = self.Cmq
        Cmq_y = self.Cmq
        op_check = 0    # open parachute flag (0:No 1:Open)

        #----------------------------------------------------------------
        if self.op_flg == 1:

            #----------------------------------------------------------------
            if self.Cd_para_1st > 0.0:
                if self.op_type_1st == 0:       # Detect trajectory peak
                    if vel[2] < 0.0 and time > (self.peak_time + self.delay_time_1st) and self.peak_time != 0:
                        op_check = 1
                        Cd_para = self.Cd_para_1st

                elif self.op_type_1st == 1:     # Open by fixed-time
                    if self.op_time_1st < time:
                        op_check = 1
                        Cd_para = self.Cd_para_1st

            #----------------------------------------------------------------
            if self.Cd_para_2nd > 0.0:
                if self.op_type_2nd == 0:       # Open by height
                    if vel[2] < 0.0 and pos[2] <= self.op_height_2nd:
                        op_check = 1
                        Cd_para = self.Cd_para_2nd

                elif self.op_type_2nd == 1:     # Open by fixed-time
                    if (self.op_time_2nd + self.peak_time) < time and self.peak_time != 0:
                        op_check = 1
                        Cd_para = self.Cd_para_2nd

            #----------------------------------------------------------------
            if self.Cd_para_3rd > 0.0:
                if self.op_type_3rd == 0:       # Open by height
                    if vel[2] < 0.0 and pos[2] <= self.op_height_3rd:
                        op_check = 1
                        Cd_para = self.Cd_para_3rd

                elif self.op_type_3rd == 1:     # Open by fixed-time
                    if (self.op_time_3rd+self.peak_time) < time  and self.peak_time != 0:
                        op_check = 1
                        Cd_para = self.Cd_para_3rd

        else:
            # Coasting
            Cd = self.Cd
            Cn_p = self.Cna * alpha
            Cn_y = self.Cna * beta
            Cmq_p = self.Cmq
            Cmq_y = self.Cmq


        # Thrust condition ------------------------------------------------
        if time < self.comb_time:
            self.force_th = self.linear_thrust(time)

        #elif self.iject_2nd < time < self.iject_2nd+self.comb_time_2nd:
        #    self.force_th = self.linear_thrust_2nd(time - self.iject_2nd)

        else:
            self.force_th = 0.0

        # Dynamics --------------------------------------------------------
        if self.op_flg == 0 or op_check == 0:

            # Define force in body coordinate
            force_b[0] = self.force_th - 0.5 * dens * (vel_norm**2) * self.ref_sq * Cd
            force_b[1] = -0.5 * dens * (vel_norm**2) * self.ref_sq * Cn_y
            force_b[2] = -0.5 * dens * (vel_norm**2) * self.ref_sq * Cn_p

            # Define damping moment in body coordinate
            mom_damp[0] = 0.0
            mom_damp[1] = 0.25 * dens * vel_norm * self.ref_sq * (self.ref_len**2) * Cmq_p * omega[1]
            mom_damp[2] = 0.25 * dens * vel_norm * self.ref_sq * (self.ref_len**2) * Cmq_y * omega[2]

            # Define moment in body coordinate
            mom_b[0] = mom_damp[0]
            mom_b[1] = mom_damp[1] + force_b[2]*dif_len
            mom_b[2] = mom_damp[2] - force_b[1]*dif_len

            # TODO: Add thrust jet damping moment

            # Joint gravity force with body force
            force_b = force_b + mass * (dcm @ Gravity)


        # Derivative value condition --------------------------------------
        if nplin.norm(pos) <= self.rail_len and vel[2] >= 0.0:

            # Launcher constraint
            if force_b[0] < 0.0:
                dVar[4:17] = 0.0

            else:
                force_b[1:3] = 0.0

                dVar[4:7] = var[7:10]
                dVar[7:10] = (dcm.T @ force_b)/mass
                dVar[10:17] = 0.0


        elif self.op_flg == 1 and op_check == 1:
            # Open 1st parachute

            # Define force in body coordinate
            vel_para = vel - wind_vel
            vel_norm_para = nplin.norm(vel_para)
            drag_para = 0.5 * dens * (vel_norm_para**2) * self.para_sq * Cd_para

            dVar[4:7] = var[7:10]
            dVar[7:10] = -1.0 * drag_para * (vel_para/vel_norm_para) + mass * Gravity
            dVar[10:17] = 0.0

        elif pos[2] < -10.0:
            # Escape odeint warning
            dVar[7:10] = 0.0

        else:
            # Coasting
            dVar[4:7] = var[7:10]
            dVar[7:10] = (dcm.T @ force_b)/mass
            dVar[10] = mom_b[0] / Ix
            dVar[11] = mom_b[1] / Iyz
            dVar[12] = mom_b[2] / Iyz
            dVar[13:17] = qt.omega2dQuat(omega, quat)

        return dVar


    #########################################################################
    def execute(self, timeVec):
        """ Use odeint solver """

        # Initial value setting
        param0 = np.array([self.mass_i, self.CGlen_i, self.Iyz_i, self.Ix_i])
        pos0 = np.zeros(3)
        vel0 = np.zeros(3)
        omega0 = np.zeros(3)
        quat0 = qt.int_quat(self.rail_elev, self.rail_azi)
        self.var0 = np.r_[param0, pos0, vel0, omega0, quat0]

        raw_result = sp.integrate.odeint(self.rocket_dynamics, self.var0, timeVec,
                                         rtol=1.e-4, atol=1.e-4)

        # Slice result data
        slice_sgn = np.sign(raw_result[:,6])
        slice_step = np.argmin(slice_sgn)
        finish_flg = np.min(slice_sgn)

        if finish_flg != -1:
            print("")
            print("Computation is not finished!")
            print("Set end_time longer.")
            print("")
            sys.exit()

        slice_result = raw_result[0:slice_step, :]

        result = np.zeros([slice_result[:,0].size, slice_result[0,:].size+1])
        result[:,0:-1] = slice_result
        result[:,-1] = timeVec[0:slice_step]

        return result
