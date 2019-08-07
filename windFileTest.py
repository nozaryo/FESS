import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

alt = 10000;
inputWind = pd.read_csv('./input/wind_data/wind_august_noshiro_nominal.csv') #xxxは適当な変数
altList = inputWind['alt']
speedList = inputWind['speed']
dirList = inputWind['dir']

altList = altList.values.tolist()
speedList = speedList.values.tolist()
dirList = dirList.values.tolist()

altList.insert(0,0)
speedList.insert(0,speedList[0])
dirList.insert(0,dirList[0])




csvIndex = 0;
for i in range(len(altList)) :
    if alt < altList[i] :
        csvIndex = i-1;
        break;
    csvIndex = i-1;

vel = ( (speedList[csvIndex + 1] - speedList[csvIndex]) / (altList[csvIndex + 1] - altList[csvIndex]) ) * (alt - altList[csvIndex]) + speedList[csvIndex]; #データ間を線形結合

dir = ( (dirList[csvIndex + 1] - dirList[csvIndex]) / (altList[csvIndex + 1] - altList[csvIndex]) ) * (alt - altList[csvIndex]) + dirList[csvIndex]; #データ間を線形結合

plt.plot(altList, speedList)
plt.show()
