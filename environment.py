import sys
import numpy as np
import plot_process as post
import matplotlib.pyplot as plt


#########################################################################
def wind_method(height, wind):
    """
    Compute wind power method

    Parameter:
        pos[2]  : Z-position(height)
        wind[0] : base wind velocity
        wind[1] : wind direction

    """
    # Set parameter of power method
    power_coeff = 4.5     # Coefficient of power method
    height_ref = 2.0     # Base height [m]
    wind_ref = wind[0]
    wind_dir = wind[1]

    # Compute wind velocity by "power law"
    if height > 0.0:
        wind_str = wind_ref * (height/height_ref) ** (1.0/power_coeff)

    else:
        wind_str = 0.0

    wind_rad = np.mod(np.deg2rad(wind_dir+180),360)

    xwind = wind_str * np.cos(wind_rad)
    ywind = wind_str * np.sin(wind_rad)
    zwind = 0.0

    wind_vel = np.array([xwind, ywind, zwind])


    return wind_vel


#########################################################################
def gas_param(height):
    """
    Compute pressure, density and temperature by height

    Parameter:
        pos[2]:Z-position(height)

    return:
        pres : pressure [Pa]
        dens : density [kg/m^3]
        temp : temperature [K]

    """

    # Set standard value
    temp0 = 10.0    # [Celsius]
    pres0 = 1.013e5    # [Pa]
    gas_R = 287.0    # [J/kgK]
    gravity = 9.81    # [m/s^2]
    base_kelvin = 273.15
    d_temp = 6.5e-3

    # Compute boundary value
    height1 = 11000.0
    height2 = 20000.0
    height3 = 32000.0

    temp0 = base_kelvin + temp0
    temp1 = temp0 - d_temp * height1
    pres1 = pres0 * (temp1/temp0)**(gravity/(d_temp*gas_R))
    pres2 = pres1 * np.exp((gravity*(height1-height2))/(gas_R*temp1))


    # Compute temperatute and pressure by height
    if(height < height1):
        d_temp = 6.5e-3

        temp = temp0 - d_temp * height
        pres = pres0 * (temp/temp0)**(gravity/(d_temp*gas_R))

    elif(height1 <= height < height2):
        d_temp = 0.0

        temp = temp1
        pres = pres1 * np.exp((gravity*(height1-height))/(gas_R*temp1))

    elif(height2 <= height < height3):
        d_temp = -1.0e-3

        temp = temp1 - d_temp * (height-height2)
        pres = pres2 * (temp/temp1)**(gravity/(d_temp*gas_R))

    else:
        print("Not define above 32000 m")
        sys.exit()

    # Compute density
    dens = pres/(gas_R*temp)

    return temp, pres, dens


if __name__ == "__main__":

    #print("Gas Parameter check")

    #count = 1000
    #max_height = height3
    #pos = np.zeros([count,3])
    #pos[:,2] = np.linspace(0.0, max_height, count)

    #pres = np.zeros(count)
    #dens = np.zeros(count)
    #temp = np.zeros(count)

    #for i in range(count):
    #    temp[i], pres[i], dens[i] = gas_param(pos[i,2])

    #post.plot_gas_dist(pos[:,2], temp, pres, dens)

    print("Wind plot")

    count = 1000
    ref_height = 2.0
    ref_wind = np.array([3.0,180.0])
    max_height = 10.0
    height = np.zeros(count)
    wind = np.zeros([3,count])

    for i in range(count):

        height[i] = (max_height/count) * i
        wind[:,i] = wind_method(height[i],ref_wind)

    plt.plot(wind[0,:],height)
    plt.xlim(0,wind[0,-1]+1.0)
    plt.xlabel("wind velocity[m/s]", fontsize=18)
    plt.ylabel("altitude[m]", fontsize=18)
    plt.show()
