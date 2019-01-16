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

class lineFunction:
    def calBy2point(self,x1,y1,x2,y2):
        self.a = (y2 -y1)/(x2 -x1)
        self.b = (x2*y1 - x1*y2) / ( x2 - x1 )

    def y(self,x):
        return self.a * x + self.b

img_origin = np.array([277,-210]) #(x,y)
safeArea_originPx = np.array([503,-435])
mag_dec_deg = 7.53   # [deg]
pixel2meter = 7.856742013157846
border_distance = 500 # 使ってない
border_point1px = np.array([189,-389])
border_point2px = np.array([458,-119])

fig, ax = plt.subplots(figsize=(11,5))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_minor_locator(MultipleLocator(100))

color_line = '#ffbf00'
# Set map image
img_map = Image.open("./map/Izu_sea_2019_2.png")
img_map = img_map.rotate(mag_dec_deg)
img_list = np.asarray(img_map)
img_width = img_map.size[0]
img_height = img_map.size[1]

safeArea_originM = safeArea_originPx
safeArea_originM = (safeArea_originPx - img_origin)*pixel2meter
border_point1m = (border_point1px - img_origin)*pixel2meter
border_point2m = (border_point2px - img_origin)*pixel2meter

origin = np.array([0,0])

#ine_O2safeArea = lineFunction()
#line_O2safeArea.calBy2point(origin[0] ,origin[1] , safeArea_originM[0] , safeArea_originM[1])
#x = np.arange(0, 5000, 1)
#y = line_O2safeArea.y(x)

border = lineFunction()
border.calBy2point(border_point1m[0],border_point1m[1],border_point2m[0],border_point2m[1])
x = np.arange(border_point1m[0], border_point2m[0], 1)
y = border.y(x)

img_left = -1.0 * img_origin[0] * pixel2meter
img_right = (img_width - img_origin[0]) * pixel2meter
img_top = -img_origin[1] * pixel2meter
img_bottom = -1.0 * (img_height + img_origin[1]) * pixel2meter
ax.plot(safeArea_originM[0],safeArea_originM[1],"yo")
ax.plot(0,0,"ro")

ax.plot(x,y,color=color_line,linestyle="dashed")


radius_safeArea = 2500

safeArea = patches.Circle(xy=safeArea_originM, radius=radius_safeArea,ec=color_line, fill=False)
ax.add_patch(safeArea)

border_length = np.sqrt( radius_safeArea**2 - (radius_safeArea - border_distance)**2 )

ax.imshow(img_list, extent=(img_left, img_right, img_bottom, img_top))
#ax.imshow(img_list)

plt.show()
