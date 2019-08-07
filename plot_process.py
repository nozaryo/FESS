import os
import math
import datetime
import json
import numpy as np
import quaternion as qt
import environment
import numpy.linalg as nplin
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import *
from PIL import Image
from mpl_toolkits.mplot3d import Axes3D

class lineFunction:
    def calBy2point(self,x1,y1,x2,y2):
        self.a = (y2 -y1)/(x2 -x1)
        self.b = (x2*y1 - x1*y2) / ( x2 - x1 )

    def y(self,x):
        return self.a * x + self.b

class JudgeInside():
    def __init__(self):
        print("Jude inside : ON")


    def set_limit_area(self, xy_range):

        # Check range area is close or not
        if np.allclose(xy_range[0,:], xy_range[-1,:]):
            print("")
            print("Range is close.")
            print("")

            self.xy_range = xy_range

        else:
            print("")
            print("Range area is not close.")
            print("Connect first point and last point automatically.")
            print("")

            point_first = xy_range[0,:]
            self.xy_range = np.vstack((xy_range, point_first))


    def set_dengerArea(self, xy_center, lim_radius = 50.0):

        self.xy_center = xy_center
        self.lim_radius = lim_radius

    def set_safeArea(self, center, radius):
        self.radius = radius
        self.center = center

    def set_border(self,line,upFlag):
        self.border_line = line
        self.upFlag = upFlag

    def judge_inside(self, check_point, place):

        # Check limit circle area is defined
        try:
            self.xy_center
            dengerArea_flag = True

        except AttributeError:
            dengerArea_flag = False

        try:
            self.xy_range
            safeArea_flag = False

        except AttributeError:
            safeArea_flag = True

        try:
            self.border_line
            border_flag = True

        except AttributeError:
            border_flag = False


        # Initialize count of line cross number
        cross_num = 0



        # Judge inside or outside by cross number
        if safeArea_flag:
            check_point_moved = check_point - self.center
            if np.linalg.norm(check_point_moved) <= self.radius:
                cross_num = 1
                if border_flag:
                    if self.upFlag:
                        if self.border_line.y(check_point[0]) > check_point[1]:
                            cross_num = 0
                    else:
                        if self.border_line.y(check_point[0]) < check_point[1]:
                            cross_num = 0
            else:
                cross_num = 0

        else:
            # Count number of range area
            point_num = self.xy_range.shape[0]

            if place == 'no_place':
                cross_num = 1
            else:
                for point in range(point_num - 1):

                    point_ymin = np.min(self.xy_range[point:point+2, 1])
                    point_ymax = np.max(self.xy_range[point:point+2, 1])

                    if check_point[1] == self.xy_range[point, 1]:

                        if check_point[0] < self.xy_range[point, 0]:
                            cross_num += 1

                        else:
                            pass

                    elif point_ymin < check_point[1] < point_ymax:

                        dx = self.xy_range[point+1, 0] - self.xy_range[point, 0]
                        dy = self.xy_range[point+1, 1] - self.xy_range[point, 1]

                        if dx == 0.0:
                            # Line is parallel to y-axis
                            judge_flag = self.xy_range[point, 1] - check_point[1]

                        elif dy == 0.0:
                            # Line is parallel to x-axis
                            judge_flag = -1.0

                        else:
                            # y = ax + b (a:slope,  b:y=intercept)
                            slope = dy / dx
                            y_intercept = self.xy_range[point, 1] - slope * self.xy_range[point, 0]

                            # left:y,  right:ax+b
                            left_eq = check_point[1]
                            right_eq = slope * check_point[0] + y_intercept

                            judge_flag = slope * (left_eq - right_eq)


                        if judge_flag > 0.0:
                            # point places left side of line
                            cross_num += 1.0

                        elif judge_flag < 0.0:
                            # point places right side of line
                            pass

                    else:
                        pass
        # odd number : inside,  even number : outside
        judge_result = np.mod(cross_num, 2)
        #print(cross_num)
        #print(judge_result)

        # check judge circle mode
        if dengerArea_flag == True:

            center_num = self.xy_center.shape[0]

            for center in range(center_num):
                # Compute distance between drop_point and center of limit circle
                length_point = np.sqrt((check_point[0] - self.xy_center[center, 0])**2 + \
                                       (check_point[1] - self.xy_center[center, 1])**2)

                # Judge in limit circle or not
                if length_point <= self.lim_radius:
                    judge_result = np.bool(False)

                else:
                    pass

        else:
            pass

        # Convert from float to bool (True:inside,  False:outside)
        judge_result = np.bool(judge_result)

        return judge_result


class PlotProcess():
    """
    Compute physical values in this method.

    input:
        result : rocket dynamics output
        wind_cond : wind condition (0:veloity, 1:direction)

    """

    def __init__(self):
        """
        Compute physical values in this method.

        input:
            result : rocket dynamics output
            wind_cond : wind condition (0:veloity, 1:direction)

        """


        # Define constant value
        self.gas_R = 287
        self.gas_gamma = 1.4

        # Initial case number
        self.case = -1

        # Generate empty matrix to store the data
        max_case = 1000

        self.wind_cond = np.zeros([max_case, 2])
        self.drop_point = np.zeros([max_case, 2])
        self.max_height = np.zeros(max_case)
        self.max_vel = np.zeros(max_case)
        self.max_mach = np.zeros(max_case)
        self.max_kakudo_deg = np.zeros(max_case)
        self.launch_clear_time = np.zeros(max_case)
        self.launch_clear_vel = np.zeros(max_case)
        self.drop_vel = np.zeros(max_case)



    def set_variety(self, result ,wind_cond,launcher_condition):
        env = environment.setEnv()
        env.wind_setting()



        #print("post processing...")

        # set launcher condition
        self.rail_len = launcher_condition[0]
        self.rail_elev = launcher_condition[1]
        self.rail_azi = launcher_condition[2]

        # Set total step number
        total_step = result.shape[0]
        self.step_num = total_step

        # Count the number of cases
        self.case += 1
        #print("Case number : {0}".format(self.case))

        # Generate matrix
        self.wind = np.zeros([total_step, 3])
        self.angle_body_rad = np.zeros([total_step, 2])
        self.dcm_body = np.zeros([3, 3, total_step])
        self.vel_body = np.zeros([total_step, 3])
        self.body_dir = np.zeros([total_step, 3])

        # Allocate value from result
        self.mass = result[:,0]
        self.len_CG = result[:,1]
        self.Iyz = result[:,2]
        self.Ix = result[:,3]
        self.pos = result[:,4:7]
        self.height = result[:,6]
        self.drop_point[self.case,:] = self.pos[-1,0:2]
        self.vel = result[:,7:10]
        self.omega = result[:,10:13]
        self.quat = result[:,13:17]
        self.time_vec = result[:,-1]
        self.dt = self.time_vec[1] - self.time_vec[0]

        # Compute vector length
        self.pos_norm = nplin.norm(self.pos, axis=1)
        self.vel_norm = nplin.norm(self.vel, axis=1)
        self.quat_norm = nplin.norm(self.quat[:,0:4], axis=1)
        self.drop_vel[self.case] = self.vel_norm[-1]

        # Compute height related parameter
        for step in range(total_step):

            # Gas parameter
            self.temp, self.pres, self.dens = env.gas_param(self.height[step])

            # Wind parameter
            self.wind[step,:] = env.wind_method(self.height[step], wind_cond)
            self.dcm_body[:,:,step] = qt.quat2dcm(self.quat[step,:])
            self.vel_body[step,:] = self.dcm_body[:,:,step] @ (self.vel[step,:] - self.wind[step,:])

            # Launch clear time/velocity
            if self.pos_norm[step] <= self.rail_len:
                self.launch_clear_time[self.case] = step * self.dt
                self.launch_clear_vel[self.case] = self.vel_norm[step]

                # Body angle constraint
                self.angle_body_rad[step,:] = 0.0

            else:
                # Compute body angle
                self.angle_body_rad[step,0] = np.arctan(self.vel_body[step,2]/self.vel_body[step,0])
                self.angle_body_rad[step,1] = np.arctan(self.vel_body[step,1]/self.vel_body[step,0])

        #
        self.sonic_speed = np.sqrt(self.gas_gamma * self.gas_R * self.temp)
        self.vel_mach = self.vel_norm / self.sonic_speed
        self.angle_body_deg = np.rad2deg(self.angle_body_rad)

        # Store wind condition
        self.wind_cond[self.case, :] = wind_cond

        theta = -1*np.arctan(self.quat[:,0]/self.quat[:,2])
        self.body_dir[:,0] = np.cos(theta)
        self.body_dir[:,1] = 0
        self.body_dir[:,2] = np.sin(theta)

        u = self.quat[:,0]
        v = self.quat[:,1]
        w = self.quat[:,2]
        sin = np.sin(self.quat[:,3])
        cos = np.cos(self.quat[:,3])
        rot_body = np.array([
        [u**2*(1-cos)+cos , u*v*(1-cos)-w*sin , w*u*(1-cos)+v*sin ],
        [u*v*(1-cos)+w*sin , v**2*(1-cos)+cos , v*w*(1-cos)-u*sin ],
        [w*u*(1-cos)-v*sin , v*w*(1-cos)+u*sin, w**2*(1-cos)+cos ]])

        for i in range(self.body_dir.shape[0]) :
            self.body_dir[i,:] = rot_body[:,:,i] @ np.array(self.body_dir[i,:])

        max_vel_index = np.argmax(self.vel_norm)
        size = (self.vel_norm[max_vel_index] * self.quat_norm[max_vel_index])

        if size == 0 :
            kakudo = 0
        else :
            kakudo = np.arccos(np.dot(np.array(self.quat[max_vel_index,1:4]),np.array(self.vel[max_vel_index,:]))/size)


        # Compute max number
        self.max_height[self.case] = np.max(self.height)
        self.max_vel[self.case] = np.max(self.vel_norm)
        self.max_mach[self.case] = np.max(self.vel_mach)
        self.max_kakudo_deg[self.case] = np.rad2deg(kakudo)


    def plot_detail(self):
        """ Plot result in detail """

        # Show results
        print("Launch clear time     : {0} sec".format(self.launch_clear_time[0]))
        print("Launch clear velocity : {0} m/s".format(self.launch_clear_vel[0]))
        print("")
        print("max height   : {0} m".format(self.max_height[0]))
        print("max velocity : {0} m/s".format(self.max_vel[0]))
        print("max Mach     : {0}".format(self.max_mach[0]))
        #print("max kakudo between axis and vel: {0} deg".format(self.max_kakudo_deg[0]))
        print("drop point   : {0}".format(self.drop_point[0,:]))

#        plt.figure()
#        plt.plot(self.time_vec, self.wind, label='vel_norm')
#        plt.xlabel("time[sec]")
#        plt.ylabel("vel[m/s]")
#        plt.legend()
#        plt.show()

#        plt.show()
##
##        fig = plt.figure()
##        ax.plot(self.pos[:,0], self.pos[:,1], self.pos[:,2])
##
##        range_lim = np.max(np.absolute(self.pos))
##        ax.set_xlim(-range_lim,range_lim)
##        ax.set_ylim(-range_lim,range_lim)
##        ax.set_zlim(0,)
##
##        ax.set_xlabel("X[m]")
#        ax.set_ylabel("Y[m]")
#        ax.set_zlabel("Up[m]")
##
##        plt.show()

        self.ax.plot(self.pos[-1,0], self.pos[-1,1],
                 color='red',marker='o')
        plt.xlabel("x(m)")
        plt.ylabel("y(m)")
        plt.legend()
        plt.grid(which='both',linestyle=":")
        self.ax.set_aspect('equal')
        now = datetime.datetime.today()
        plt.savefig('./result/detail_{}_{}_{}_{}_{}_{}.png'.format(now.year,now.month,now.day,now.hour,now.minute,now.second))
        plt.show();

        with open("./result/result.csv", "w") as res:
            res.write(str(datetime.datetime.today())+'\n')
            res.write("time[sec],mass[kg],len_CG[m],Iyz[kg*m^2],Ix[kg*m^2],pos_x[m],pos_y[m],pos_z[m],vel_x[m/s],vel_y[m/s],vel_z[m/s],roll,pitch,yaw,quat_x,quat_y,quat_z,quat_w\n")
            for step in range(self.step_num):
                res.write(str(self.time_vec[step])+',')
                res.write(str(self.mass[step])+',')
                res.write(str(self.len_CG[step])+',')
                res.write(str(self.Iyz[step])+',')
                res.write(str(self.Ix[step])+',')
                res.write(str(self.pos[step,0])+',')
                res.write(str(self.pos[step,1])+',')
                res.write(str(self.pos[step,2])+',')
                res.write(str(self.vel[step,0])+',')
                res.write(str(self.vel[step,1])+',')
                res.write(str(self.vel[step,2])+',')
                res.write(str(self.omega[step,0])+',')
                res.write(str(self.omega[step,1])+',')
                res.write(str(self.omega[step,2])+',')
                res.write(str(self.quat[step,0])+',')
                res.write(str(self.quat[step,1])+',')
                res.write(str(self.quat[step,2])+',')
                res.write(str(self.quat[step,3])+'\n')



    def set_map(self, place,rail_cond):
        """ Set lat/lon coordinates to define MAP """

        earth_radius = 6378150.0    # [km]

        self.place = place

        if place == 'Izu_land':

            self.fig, self.ax = plt.subplots(figsize=(15,8))
            self.ax.xaxis.set_minor_locator(MultipleLocator(100))
            self.ax.yaxis.set_minor_locator(MultipleLocator(100))

            # Set lat/long coordinates
            # point_origin : map origin
            # point_center : circle limt area
            # point_range  : limit area vertex
            self.lim_radius = 50.0   # define circle limit area

            self.point_origin = np.array([34.735972, 139.420944])

            self.point_center = np.array([[34.735972, 139.420944],
                                          [34.735390, 139.421377],
                                          [34.731230, 139.423150]])

            self.point_range = np.array([[34.735715, 139.420922],
                                         [34.731750, 139.421719],
                                         [34.733287, 139.424590],
                                         [34.736955, 139.426038],
                                         [34.738908, 139.423597],
                                         [34.740638, 139.420681],
                                         [34.741672, 139.417387],
                                         [34.735715, 139.420922],
                                         ])

            self.point_center_rel = self.point_center - self.point_origin
            self.point_range_rel = self.point_range - self.point_origin

            # Set magnetic declination
            self.mag_dec_deg = -7.53   # [deg]

            mag_dec_rad = np.deg2rad(self.mag_dec_deg)
            mat_rot = np.array([[np.cos(mag_dec_rad), -1 * np.sin(mag_dec_rad)],
                                [np.sin(mag_dec_rad), np.cos(mag_dec_rad)]])

            # Convert lat/lon to meter
            self.lat2met = 2 * math.pi * earth_radius / 360.0
            self.lon2met = 2 * math.pi * earth_radius * np.cos(np.deg2rad(self.point_origin[0])) / 360.0

            # Convert from lat/long to meter (ENU coordinate)
            self.xy_center = np.zeros(self.point_center.shape)
            self.xy_range = np.zeros(self.point_range.shape)

            self.xy_center[:,0] = self.lon2met * self.point_center_rel[:,1]
            self.xy_center[:,1] = self.lat2met * self.point_center_rel[:,0]
            self.xy_range[:,0] = self.lon2met * self.point_range_rel[:,1]
            self.xy_range[:,1] = self.lat2met * self.point_range_rel[:,0]

            # Apply magnetic effect
            for i in range(self.point_center.shape[0]):
                self.xy_center[i,:] = mat_rot @ self.xy_center[i,:]

            for i in range(self.point_range.shape[0]):
                self.xy_range[i,:] = mat_rot @ self.xy_range[i,:]

            # Setup MAP image --------------------------
            # Convert pixel to meter
            pixel2meter = 0.946981208125

            # Set map image
            img_map = Image.open("./map/Izu_map_mag.png")
            img_list = np.asarray(img_map)
            img_height = img_map.size[0]
            img_width = img_map.size[1]
            img_origin = np.array([722, 749])    # TODO : compute by lat/long of launcher point

            # Define image range
            img_left =   -1.0 * img_origin[0] * pixel2meter
            img_right = (img_width - img_origin[0]) * pixel2meter
            img_top = img_origin[1] * pixel2meter
            img_bottom = -1.0 * (img_height - img_origin[1]) * pixel2meter

            #plot.figure(figsize=(10,8))
            self.ax.imshow(img_list, extent=(img_left, img_right, img_bottom, img_top))

            # Define color
            color_line = '#ffbf00'    # Yellow
            color_circle = 'r'    # Red

            # Set circle object
            #self.ax = plt.axes()

            # plot limit area
            for i in range(self.point_center.shape[0]):
                circle = patches.Circle(xy=self.xy_center[i,:], radius=self.lim_radius,
                                        ec=color_circle, fill=False)
                self.ax.add_patch(circle)
                self.ax.plot(self.xy_center[i,0], self.xy_center[i,1], '.', color=color_circle)


            self.ax.plot(self.xy_range[:,0], self.xy_range[:,1], '--', color=color_line)
            #plt.plot(self.xy_range[:,0], self.xy_range[:,1], '.', color='r')



        elif place == 'nosiro_land':

            # Set lat/long coordinates
            # point_origin : map origin
            # point_center : circle limt area
            # point_range  : limit area vertex
            self.fig, self.ax = plt.subplots(figsize=(10,8))
            self.ax.xaxis.set_minor_locator(MultipleLocator(100))
            self.ax.yaxis.set_minor_locator(MultipleLocator(100))

            self.lim_radius = 5.0   # define circle limit area

            self.point_origin = np.array([40.138633, 139.984850])

            self.point_center = np.array([[40.138633, 139.984850]])

            self.point_range = np.array([[40.139725, 139.983939],
                                         [40.136127, 139.982133],
                                         [40.135607, 139.981753],
                                         [40.134911, 139.981451],
                                         [40.134821, 139.981692],
                                         [40.135639, 139.983324],
                                         [40.137052, 139.984608],
                                         [40.138053, 139.985781],
                                         [40.139075, 139.986297],
                                         [40.139725, 139.983939]])

            self.point_center_rel = self.point_center - self.point_origin
            self.point_range_rel = self.point_range - self.point_origin

            # Set magnetic declination
            self.mag_dec_deg = -8.9   # [deg]

            mag_dec_rad = np.deg2rad(self.mag_dec_deg)
            mat_rot = np.array([[np.cos(mag_dec_rad), -1 * np.sin(mag_dec_rad)],
                                [np.sin(mag_dec_rad), np.cos(mag_dec_rad)]])

            # Convert lat/lon to meter
            self.lat2met = 2 * math.pi * earth_radius / 360.0
            self.lon2met = 2 * math.pi * earth_radius * np.cos(np.deg2rad(self.point_origin[0])) / 360.0

            # Convert from lat/long to meter (ENU coordinate)
            self.xy_center = np.zeros(self.point_center.shape)
            self.xy_range = np.zeros(self.point_range.shape)

            self.xy_center[:,0] = self.lon2met * self.point_center_rel[:,1]
            self.xy_center[:,1] = self.lat2met * self.point_center_rel[:,0]
            self.xy_range[:,0] = self.lon2met * self.point_range_rel[:,1]
            self.xy_range[:,1] = self.lat2met * self.point_range_rel[:,0]

            self.xy_range[0,]=[-83+12,113+7]
            self.xy_range[1,]=[-254+12,-250+7]
            self.xy_range[2,]=[-288.0+12,-304.0+7]
            self.xy_range[3,]=[-320.0+12,-380.0+7]
            self.xy_range[4,]=[-296.0+12,-388.0+7]
            self.xy_range[5,]=[-146.0+12,-300.0+7]
            self.xy_range[6,]=[-29+12,-155+7]
            self.xy_range[7,]=[72+12,-52+7]
            self.xy_range[8,]=[125+12,52+7]
            self.xy_range[9,]=self.xy_range[0,]

            # Apply magnetic effect
            for i in range(self.point_center.shape[0]):
                self.xy_center[i,:] = mat_rot @ self.xy_center[i,:]

            for i in range(self.point_range.shape[0]):
                self.xy_range[i,:] = mat_rot @ self.xy_range[i,:]


            # Setup MAP image --------------------------
            # Convert pixel to meter
            pixel2meter = 0.58479532

            # Set map image
            img_map = Image.open("./map/nosiro2.png")
            img_map = img_map.rotate(self.mag_dec_deg)
            img_list = np.asarray(img_map)
            img_height = img_map.size[0]
            img_width = img_map.size[1]
            img_origin = np.array([690-82,1293+173])  # TODO : compute by lat/long of launcher point

            # Define image range
            img_left =   -1.0 * img_origin[1] * pixel2meter
            img_right = (img_width - img_origin[1]) * pixel2meter
            img_top = img_origin[0] * pixel2meter
            img_bottom = -1.0 * (img_height - img_origin[0]) * pixel2meter
            img_top = img_top*2.0
            img_bottom = img_bottom*2.0
            img_left = img_left/2.0
            img_right = img_right/2.0
            #plot.figure(figsize=(10,8))
            self.ax.imshow(img_list, extent=(img_left, img_right, img_bottom, img_top))

            # Define color
            color_line = '#ffbf00'    # Yellow
            color_circle = 'r'    # Red

            # Set circle object
            #self.ax = plt.axes()

            # plot limit area
            for i in range(self.point_center.shape[0]):
                circle = patches.Circle(xy=self.xy_center[i,:], radius=self.lim_radius,
                                        ec=color_circle, fill=False)
                self.ax.add_patch(circle)
                self.ax.plot(self.xy_center[i,0], self.xy_center[i,1], '.', color=color_circle)


            self.ax.plot(self.xy_range[:,0], self.xy_range[:,1], '--', color=color_line)
            #plt.plot(self.xy_range[:,0], self.xy_range[:,1], '.', color='r')

        elif place == 'no_place':

            # Set lat/long coordinates
            # point_origin : map origin
            # point_center : circle limt area
            # point_range  : limit area vertex
            self.fig, self.ax = plt.subplots(figsize=(11,5))
            self.ax.xaxis.set_minor_locator(MultipleLocator(100))
            self.ax.yaxis.set_minor_locator(MultipleLocator(100))
            self.lim_radius = 10.0   # define circle limit area


            # Set magnetic declination

            # Convert from lat/long to meter (ENU coordinate)
            self.xy_center = np.zeros((1,2))
            self.xy_range = np.zeros((5,2))

            #self.xy_range[0,]=[100,100]
            #self.xy_range[1,]=[100,-100]
            #self.xy_range[2,]=[-100,-100]
            #self.xy_range[3,]=[-100,100]
            #self.xy_range[4,]=self.xy_range[0,]

            # Convert pixel to meter
            pixel2meter = 0.34862385

            # Set map image
            img_map = Image.open("./map/baseboll.png")
            img_map = img_map.rotate(rail_cond[2]-90.0)
            #img_map = img_map.rotate(self.mag_dec_deg)
            img_list = np.asarray(img_map)

            img_height = img_map.size[1]
            img_width = img_map.size[0]
            origin = np.array([250-img_width/2.0 ,380-img_height/2.0])  # TODO : compute by lat/long of launcher point
            origin = np.array([origin[0]*np.cos(np.deg2rad(-(rail_cond[2]-90.0)))-origin[1]*np.sin(np.deg2rad(-(rail_cond[2]-90.0))),
                                   origin[0]*np.sin(np.deg2rad(-(rail_cond[2]-90.0)))+origin[1]*np.cos(np.deg2rad(-(rail_cond[2]-90.0)))])
            img_origin = np.array([origin[0]+img_width/2.0,origin[1]+img_height/2.0])
            # Define image range
            img_left =   -1.0 * img_origin[0] * pixel2meter
            img_right = (img_width - img_origin[0]) * pixel2meter
            img_top = img_origin[1] * pixel2meter
            img_bottom = -1.0 * (img_height - img_origin[1]) * pixel2meter

            self.ax.imshow(img_list, extent=(img_left, img_right, img_bottom, img_top))
            # Setup MAP image --------------------------
            # Define color
            color_line = '#ffbf00'    # Yellow
            color_circle = 'r'    # Red

            # Set circle object
            #self.ax = plt.axes()
            circle = patches.Circle(xy=[0,0], radius=self.lim_radius,
                                        ec=color_circle, fill=False)
            self.ax.add_patch(circle)
            self.ax.plot(0, 0, '.', color=color_circle)

            #self.ax.plot(self.xy_range[:,0], self.xy_range[:,1], '--', color=color_line)
            #plt.plot(self.xy_range[:,0], self.xy_range[:,1], '.', color='r')

        elif place == 'nosiro_land2': #何で作ったか思い出せない

            # Set lat/long coordinates
            # point_origin : map origin
            # point_center : circle limt area
            # point_range  : limit area vertex
            self.fig, self.ax = plt.subplots(figsize=(10,8))
            self.ax.xaxis.set_minor_locator(MultipleLocator(100))
            self.ax.yaxis.set_minor_locator(MultipleLocator(100))

            self.lim_radius = 5.0   # define circle limit area

            self.point_origin = np.array([40.138633, 139.984850])

            self.point_center = np.array([[40.138633, 139.984850]])

            self.point_range = np.array([[40.139725, 139.983939],
                                         [40.136127, 139.982133],
                                         [40.135607, 139.981753],
                                         [40.134911, 139.981451],
                                         [40.134821, 139.981692],
                                         [40.135639, 139.983324],
                                         [40.137052, 139.984608],
                                         [40.138053, 139.985781],
                                         [40.139075, 139.986297],
                                         [40.139725, 139.983939]])

            self.point_center_rel = self.point_center - self.point_origin
            self.point_range_rel = self.point_range - self.point_origin

            # Set magnetic declination
            self.mag_dec_deg = -7.53   # [deg]

            mag_dec_rad = np.deg2rad(self.mag_dec_deg)
            mat_rot = np.array([[np.cos(mag_dec_rad), -1 * np.sin(mag_dec_rad)],
                                [np.sin(mag_dec_rad), np.cos(mag_dec_rad)]])

            # Convert lat/lon to meter
            self.lat2met = 2 * math.pi * earth_radius / 360.0
            self.lon2met = 2 * math.pi * earth_radius * np.cos(np.deg2rad(self.point_origin[0])) / 360.0

            # Convert from lat/long to meter (ENU coordinate)
            self.xy_center = np.zeros(self.point_center.shape)
            self.xy_range = np.zeros(self.point_range.shape)

            self.xy_center[:,0] = self.lon2met * self.point_center_rel[:,1]
            self.xy_center[:,1] = self.lat2met * self.point_center_rel[:,0]
            self.xy_range[:,0] = self.lon2met * self.point_range_rel[:,1]
            self.xy_range[:,1] = self.lat2met * self.point_range_rel[:,0]

            # Apply magnetic effect
            for i in range(self.point_center.shape[0]):
                self.xy_center[i,:] = mat_rot @ self.xy_center[i,:]

            for i in range(self.point_range.shape[0]):
                self.xy_range[i,:] = mat_rot @ self.xy_range[i,:]

            # Setup MAP image --------------------------


            # Setup MAP image --------------------------
            # Convert pixel to meter
            pixel2meter = 0.58479532

            # Set map image
            img_map = Image.open("./map/nosiro2.png")
            img_map = img_map.rotate(self.mag_dec_deg)
            img_list = np.asarray(img_map)
            img_height = img_map.size[0]
            img_width = img_map.size[1]
            img_origin = np.array([690,1293])  # TODO : compute by lat/long of launcher point

            # Define image range
            img_left =   -1.0 * img_origin[1] * pixel2meter
            img_right = (img_width - img_origin[1]) * pixel2meter
            img_top = img_origin[0] * pixel2meter
            img_bottom = -1.0 * (img_height - img_origin[0]) * pixel2meter

            #plot.figure(figsize=(10,8))
            self.ax.imshow(img_list, extent=(img_left, img_right, img_bottom, img_top))

            # Define color
            color_line = '#ffbf00'    # Yellow
            color_circle = 'r'    # Red

            # Set circle object
            #self.ax = plt.axes()

            # plot limit area
            for i in range(self.point_center.shape[0]):
                circle = patches.Circle(xy=self.xy_center[i,:], radius=self.lim_radius,
                                        ec=color_circle, fill=False)
                self.ax.add_patch(circle)
                self.ax.plot(self.xy_center[i,0], self.xy_center[i,1], '.', color=color_circle)


            self.ax.plot(self.xy_range[:,0], self.xy_range[:,1], '--', color=color_line)
            #plt.plot(self.xy_range[:,0], self.xy_range[:,1], '.', color='r')

        elif place == "Izu_sea":

            self.fig, self.ax = plt.subplots(figsize=(13,8))
            self.ax.xaxis.set_minor_locator(MultipleLocator(500))
            self.ax.yaxis.set_minor_locator(MultipleLocator(500))

            img_origin = np.array([360,-159]) #(x,y)
            safeArea_originPx = np.array([520,-434])
            self.mag_dec_deg = -7.53   # [deg]
            pixel2meter = 7.856742013157846
            border_distance = 500 # 使ってない
            border_point1px = np.array([227,-309])
            border_point2px = np.array([557,-117])

            color_line = '#ffbf00'
            # Set map image
            img_map = Image.open("./map/Izu_sea_2019_2.png")
            img_map = img_map.rotate(self.mag_dec_deg)
            img_list = np.asarray(img_map)
            img_width = img_map.size[0]
            img_height = img_map.size[1]

            self.safeArea_originM = (safeArea_originPx - img_origin)*pixel2meter
            border_point1m = (border_point1px - img_origin)*pixel2meter
            border_point2m = (border_point2px - img_origin)*pixel2meter

            origin = np.array([0,0])
            self.xy_center = np.zeros([1,2])

            #ine_O2safeArea = lineFunction()
            #line_O2safeArea.calBy2point(origin[0] ,origin[1] , self.safeArea_originM[0] , self.safeArea_originM[1])
            #x = np.arange(0, 5000, 1)
            #y = line_O2safeArea.y(x)

            self.border = lineFunction()
            self.border.calBy2point(border_point1m[0],border_point1m[1],border_point2m[0],border_point2m[1])
            x = np.arange(border_point1m[0], border_point2m[0], 1)
            y = self.border.y(x)

            img_left = -1.0 * img_origin[0] * pixel2meter
            img_right = (img_width - img_origin[0]) * pixel2meter
            img_top = -img_origin[1] * pixel2meter
            img_bottom = -1.0 * (img_height + img_origin[1]) * pixel2meter
            self.ax.plot(self.safeArea_originM[0],self.safeArea_originM[1],"yo")
            self.ax.plot(0,0,"ro")

            self.ax.plot(x,y,color=color_line,linestyle="dashed")


            self.radius_safeArea = 2500

            safeArea = patches.Circle(xy=self.safeArea_originM, radius=self.radius_safeArea,ec=color_line, fill=False)
            self.ax.add_patch(safeArea)


            self.ax.imshow(img_list, extent=(img_left, img_right, img_bottom, img_top))
            #ax.imshow(img_list)

        elif place == 'nosiro_sea':

            # Set lat/long coordinates
            # point_origin : map origin
            # point_center : circle limt area
            # point_range  : limit area vertex
            self.fig, self.ax = plt.subplots(figsize=(11,5))
            self.ax.xaxis.set_minor_locator(MultipleLocator(500))
            self.ax.yaxis.set_minor_locator(MultipleLocator(500))
            self.lim_radius = 10.0   # define circle limit area
            self.mag_dec_deg = -8.94


            # Set magnetic declination

            # Convert from lat/long to meter (ENU coordinate)
            self.xy_center = np.zeros((1,2))
            self.xy_range = np.zeros((5,2))

            #self.xy_range[0,]=[100,100]
            #self.xy_range[1,]=[100,-100]
            #self.xy_range[2,]=[-100,-100]
            #self.xy_range[3,]=[-100,100]
            #self.xy_range[4,]=self.xy_range[0,]

            # Convert pixel to meter
            pixel2meter = 7.246376811

            # Set map image
            img_map = Image.open("./map/nosiro_sea_2019.png")
            img_map = img_map.rotate(self.mag_dec_deg)
            img_list = np.asarray(img_map)

            img_height = img_map.size[1]
            img_width = img_map.size[0]
            img_origin = np.array([1207,-982])
            # Define image range
            img_left = -1.0 * img_origin[0] * pixel2meter
            img_right = (img_width - img_origin[0]) * pixel2meter
            img_top = -img_origin[1] * pixel2meter
            img_bottom = -1.0 * (img_height + img_origin[1]) * pixel2meter

            self.ax.imshow(img_list, extent=(img_left, img_right, img_bottom, img_top))
            # Setup MAP image --------------------------
            # Define color
            color_line = '#ffbf00'    # Yellow
            color_circle = 'r'    # Red

            # Set circle object
            #self.ax = plt.axes()
            circle = patches.Circle(xy=[0,0], radius=self.lim_radius,
                                        ec=color_circle, fill=False)
            self.ax.add_patch(circle)
            self.ax.plot(0, 0, '.', color=color_circle)

            #self.ax.plot(self.xy_range[:,0], self.xy_range[:,1], '--', color=color_line)
            #plt.plot(self.xy_range[:,0], self.xy_range[:,1], '.', color='r')


    def plot_scatter(self, filename, wind_case,deg,ie,op_flg,elev_mode,place):

        f = open('input/'+filename, 'r')
        stdin_json = json.load(f)
        f.close()
        stdin_info = stdin_json["info"]
        stdin_rocket = stdin_json["rocket"]
        stdin_motor = stdin_json["motor"]
        stdin_env = stdin_json["environment"]
        env.wind_file_set(str(stdin_env.get("wind_file", 0)))

        env = environment.setEnv()
        env.wind_setting()

        # Define computation pattern
        vel_pat = int(wind_case[2])
        dir_pat = int(wind_case[3])
        self.judge_result = np.zeros([dir_pat, vel_pat], dtype=bool)
        judge = JudgeInside()
        # Call JudgeInside class]
        if place != "Izu_sea":
            judge.set_limit_area(self.xy_range)
            judge.set_dengerArea(self.xy_center, lim_radius=self.lim_radius)
        else:
            judge.set_safeArea(self.safeArea_originM,self.radius_safeArea)
            judge.set_border(self.border,False)
        # Judge landing point is inside limit area or not
        for dir in range(dir_pat):
            for vel in range(vel_pat):
                case = vel * dir_pat + dir + ie*dir_pat*vel_pat
                self.judge_result[dir, vel] = judge.judge_inside(self.drop_point[case,:],self.place)

        # Plot scatter graph on MAP
        cmap = plt.get_cmap("spring")

        for iter in range(vel_pat):
            case_st = vel_pat*dir_pat*ie + iter * dir_pat
            case_en = case_st + dir_pat

            wind_vel = self.wind_cond[case_st, 0]
            label_name = str(wind_vel) + "m/s"

            x_point = np.r_[self.drop_point[case_st:case_en, 0], self.drop_point[case_st, 0]]
            y_point = np.r_[self.drop_point[case_st:case_en, 1], self.drop_point[case_st, 1]]
            self.ax.plot(x_point, y_point,
                     color=cmap(iter/vel_pat), label = label_name, marker='o')


        # show result
        max_height_buff = self.max_height
        max_vel_buff = self.max_vel
        launch_clear_vel_buff = self.launch_clear_vel
        max_kakudo_deg_buff = self.max_kakudo_deg
        max_drop_vel_buff = self.drop_vel

        print(self.judge_result)
        while 1 :
            if np.amax(max_height_buff) == 0 :
                break
            max_height_vel = (int)(np.argmax(max_height_buff)/wind_case[3])#*wind_case[1]+wind_case[0]
            max_height_dir = (np.argmax(max_height_buff)%wind_case[3])#*(360/wind_case[3])

            if self.judge_result[int(max_height_dir),int(int(max_height_vel)%(wind_case[2]))] == False :
                max_height_buff[np.argmax(max_height_buff)] = 0
            elif max_height_dir == 1 :
                max_height_buff[np.argmax(max_height_buff)] = 0
            else :
                max_height_vel = int(max_height_vel)%(wind_case[2])*wind_case[1]+wind_case[0]
                max_height_dir = max_height_dir*(360/wind_case[3])
                break

        print('Max height {}m (wind {}m/s {} deg )'.format(np.amax(max_height_buff),max_height_vel,max_height_dir))

        while 1 :
            if np.amax(max_vel_buff) == 0 :
                break
            max_vel_vel = (int)(np.argmax(max_vel_buff)/wind_case[3])#*wind_case[1]+wind_case[0]
            max_vel_dir = (np.argmax(max_vel_buff)%wind_case[3])#*(360/wind_case[3])

            if self.judge_result[int(max_vel_dir),int(int(max_vel_vel)%(wind_case[2]))] == False :
                max_vel_buff[np.argmax(max_vel_buff)] = 0
            elif max_vel_dir == 1 :
                max_vel_buff[np.argmax(max_vel_buff)] = 0
            else :
                max_vel_vel = int(max_vel_vel)%(wind_case[2])*wind_case[1]+wind_case[0]
                max_vel_dir = max_vel_dir*(360/wind_case[3])
                break

        print('Max veloity {}m/s (wind {}m/s {} deg)'.format(np.amax(max_vel_buff),max_vel_vel,max_vel_dir))

        while 1 :
            if np.amin(launch_clear_vel_buff[np.nonzero(launch_clear_vel_buff)]) == 999 :
                break
            min_launch_vel = (int)(np.argmin(launch_clear_vel_buff[np.nonzero(launch_clear_vel_buff)])/wind_case[3])#*wind_case[1]+wind_case[0]
            min_launch_dir = (np.argmin(launch_clear_vel_buff[np.nonzero(launch_clear_vel_buff)])%wind_case[3])#*(360/wind_case[3])

            if self.judge_result[int(min_launch_dir),int(int(min_launch_vel)%(wind_case[2]))] == False :
                launch_clear_vel_buff[np.argmin(launch_clear_vel_buff[np.nonzero(launch_clear_vel_buff)])] = 999
            elif min_launch_dir == 1 :
                launch_clear_vel_buff[np.argmin(launch_clear_vel_buff[np.nonzero(launch_clear_vel_buff)])] = 999
            else :
                min_launch_vel = int(min_launch_vel)%(wind_case[2])*wind_case[1]+wind_case[0]
                min_launch_dir = min_launch_dir*(360/wind_case[3])
                break

        print('Minimum launch clear vel {}m/s (wind {}m/s {} deg)'.format(np.amin(launch_clear_vel_buff[np.nonzero(launch_clear_vel_buff)]),min_launch_vel,min_launch_dir))

        while 1 :
            if np.amax(max_drop_vel_buff) == 0 :
                break
            max_drop_vel_vel = (int)(np.argmax(max_drop_vel_buff)/wind_case[3])#*wind_case[1]+wind_case[0]
            max_drop_vel_dir = (np.argmax(max_drop_vel_buff)%wind_case[3])#*(360/wind_case[3])

            if self.judge_result[int(max_drop_vel_dir),int(int(max_drop_vel_vel)%(wind_case[2]))] == False :
                max_drop_vel_buff[np.argmax(max_drop_vel_buff)] = 0
            elif max_drop_vel_dir == 1 :
                max_drop_vel_buff[np.argmax(max_drop_vel_buff)] = 0
            else :
                max_drop_vel_vel = int(max_drop_vel_vel)%(wind_case[2])*wind_case[1]+wind_case[0]
                max_drop_vel_dir = max_drop_vel_dir*(360/wind_case[3])
                break

        print('Max drop veloity {}m/s (wind {}m/s {} deg)'.format(np.amax(max_drop_vel_buff),max_drop_vel_vel,max_drop_vel_dir))



        if op_flg==1:
            plt.title("Parachute "+ str(deg)+"deg")
        else:
            plt.title("Trajectory "+ str(deg)+"deg")


        plt.xlabel("x(m)")
        plt.ylabel("y(m)")
        plt.legend()
        plt.grid(which='both',linestyle=":")
        self.ax.set_aspect('equal')


        azi = stdin_env.get("rail_azi")

        self.ax.arrow(0,0,600*np.cos(np.deg2rad(azi)),600*np.sin(np.deg2rad(azi)), width=200, head_width=500, head_length=300, fc='k', ec='k',alpha = 0.5)

        # fig内でのaxes座標を取得，戻り値はBbox
        ax_pos = self.ax.get_position()

        if "Cmq" not in stdin_rocket:
            Cmq = -(stdin_rocket["Cna"]/2)*( ( (stdin_rocket["CPlen"] - stdin_rocket["CGlen_f"])/ stdin_rocket["ref_len"] ) )**2
        else:
            Cmq = stdin_rocket.get("Cmq", -2.0)

        # fig内座標でテキストを表示 Bboxは Bbox.x0, Bbox.x1, Bbox.y0, Bbox.y1で座標を取得できる
        self.fig.text(ax_pos.x0 - 0.18, ax_pos.y1-0.4, "Team    : " + stdin_info.get("TEAM")+"\n"
                                                  "Name    : " + stdin_info.get("NAME")+"\n"
                                                  "Date    : " + str(datetime.date.today())+"\n\n\n"
                                                  "ref_len : " + str(stdin_rocket["ref_len"])+"[m]\n"
                                                  "diam    : " + str(stdin_rocket["diam"])+"[m]\n"
                                                  "CGlen_i : " + str(stdin_rocket["CGlen_i"])+"[m]\n"
                                                  "CGlen_f : " + str(stdin_rocket["CGlen_f"])+"[m]\n"
                                                  "mass_i  : " + str(stdin_rocket["mass_i"])+"[kg]\n"
                                                  "mass_f  : " + str(stdin_rocket["mass_f"])+"[kg]\n"
                                                  "Iyz_i   : " + str(stdin_rocket["Iyz_i"])+"[kg*m^2]\n"
                                                  "Iyz_f   : " + str(stdin_rocket["Iyz_f"])+"[kg*m^2]\n"
                                                  "CPlen   : " + str(stdin_rocket["CPlen"])+"[m]\n"
                                                  "Cd      : " + str(stdin_rocket["Cd"])+"\n"
                                                  "Cna     : " + str(stdin_rocket["Cna"])+"\n"
                                                  "Cmq     : " + "{:.4f}".format(Cmq)
                                                  )
        self.fig.text(ax_pos.x1+0.01 , ax_pos.y0,'Max height {:.3f}m (wind {}m/s {} deg )'.format(np.amax(max_height_buff),max_height_vel,max_height_dir)+"\n"
                                                             'Max veloity {:.3f}m/s (wind {}m/s {} deg)'.format(np.amax(max_vel_buff),max_vel_vel,max_vel_dir)+"\n"
                                                             'Minimum launch clear vel {:.3f}m/s (wind {}m/s {} deg)'.format(np.amin(launch_clear_vel_buff[np.nonzero(launch_clear_vel_buff)]),min_launch_vel,min_launch_dir)+"\n"
                                                             'Max drop veloity {:.3f}m/s (wind {}m/s {} deg)'.format(np.amax(max_drop_vel_buff),max_drop_vel_vel,max_drop_vel_dir),fontsize=7)

        now = datetime.datetime.today()
        plt.savefig('./result/{}deg_{}_{}_{}_{}_{}_{}.png'.format(int(deg),now.year,now.month,now.day,now.hour,now.minute,now.second))
        if elev_mode == 0:
            plt.show()


    def plot_orbit(var):
        """
        Plot the result

        """

        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot(var[:,4],var[:,5],var[:,6])

        range_lim = np.max(np.absolute(var[:,4:6]))
        ax.set_xlim(-range_lim,range_lim)
        ax.set_ylim(-range_lim,range_lim)
        ax.set_zlim(0,)

        ax.set_xlabel("X[m]")
        ax.set_ylabel("Y[m]")
        ax.set_zlabel("Up[m]")

    def cov_position(self, origin_earth_pos1,origin_earth_pos2,posX,posY):
        earth_radius = 6378150.0    # [m]
        mag_dec_rad = np.deg2rad(self.mag_dec_deg)
        mat_rot = np.array([[np.cos(mag_dec_rad), -1 * np.sin(mag_dec_rad)],
                            [np.sin(mag_dec_rad), np.cos(mag_dec_rad)]])
        pos=np.zeros(2)
        pos=np.array([posX,posY])
        pos = pos @ mat_rot
        def_pos2 = ( pos[0] * 360.0 )/ (2*np.pi*earth_radius*np.sin( np.deg2rad(90-origin_earth_pos1) ))
        def_pos1 = ( pos[1] * 360.0 )/ (2*np.pi*earth_radius)
        return (origin_earth_pos1 + def_pos1 , origin_earth_pos2 + def_pos2,def_pos1,def_pos2)




if __name__ == "__main__":

##    place = 'Izu_land'
##
##    plot_map = PostProcess()
##
##    plt.show()

    post = PlotProcess()
    post.set_map("Izu_sea",(0,0,0))
    print(post.cov_position(34.679730,139.438373,852.38985154 ,-439.864895))
