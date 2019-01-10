import os
import math
import datetime
import json
import numpy as np
import quaternion as qt
import environment as env
import numpy.linalg as nplin
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import *
from PIL import Image
from mpl_toolkits.mplot3d import Axes3D

fig, ax = plt.subplots(figsize=(11,5))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_minor_locator(MultipleLocator(100))
mag_dec_deg = 7.53   # [deg]
pixel2meter = 5.882352941176471
color_line = '#ffbf00'
# Set map image
img_map = Image.open("./map/Izu_sea_2019.png")
img_map = img_map.rotate(mag_dec_deg)
img_list = np.asarray(img_map)
img_width = img_map.size[0]
img_height = img_map.size[1]
img_origin = np.array([361,-178]) #(x,y)
px_limit_circle_origin = np.array([660,-486])
m_limit_circle_origin = px_limit_circle_origin
m_limit_circle_origin[0] = (px_limit_circle_origin[0]-img_origin[0])*pixel2meter
m_limit_circle_origin[1] = (px_limit_circle_origin[1]-img_origin[1])*pixel2meter

img_left =   -1.0 * img_origin[0] * pixel2meter
img_right = (img_width - img_origin[0]) * pixel2meter
img_top = -img_origin[1] * pixel2meter
img_bottom = -1.0 * (img_height + img_origin[1]) * pixel2meter
ax.plot(m_limit_circle_origin[0],m_limit_circle_origin[1],"yo")
ax.plot(0,0,"ro")

circle = patches.Circle(xy=m_limit_circle_origin, radius=2500,ec=color_line, fill=False)
ax.add_patch(circle)

ax.imshow(img_list, extent=(img_left, img_right, img_bottom, img_top))
#ax.imshow(img_list)

plt.show()
