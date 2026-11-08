"""
Project: Atmospheric Divergence and Vertical Velocity Analysis
Description:
This project implements atmospheric wind speed and direction analysis using Python as part of academic coursework.
"""

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

#上傳資料
uploaded = files.upload()

lev = 5
lon = 49
lat = 25

dy=6378000*1.875*math.pi/180
#計算φ(phi)
phi = np.linspace(15,60,lat)*math.pi/180
#計算Ω(omega)
omega = 2*math.pi/(86400)
# 計算科氏力
f = []
for i in range(25):
    f.append(2*omega*math.sin(phi[i]))

# 差分法
# x方向上的差分
def ddx(ver,mlev,mlat,mlon):
    varx = 0
    dx = dy * math.cos(phi[mlat])
    if (mlon == 0):
        # 前邊界，故使用向前差分法
        varx = (ver[mlev,mlat,mlon+1] - ver[mlev,mlat,mlon]) / dx
        return varx
    elif (mlon == lon-1):
        # 後邊界，故使用向後差分法
        varx = (ver[mlev,mlat,mlon] - ver[mlev,mlat,mlon-1]) / dx
        return varx
    else:
        # 其他非邊界項，使用中差法
        varx = (ver[mlev,mlat,mlon+1] - ver[mlev,mlat,mlon-1]) / (2*dx)
        return varx

#y方向上的差分
def ddy(ver,mlev,mlat,mlon):
    vary = 0
    if (mlat == 0):
        # 前邊界，故使用向前差分法
        vary= (ver[mlev,mlat+1,mlon] - ver[mlev,mlat,mlon]) / dy
        return vary
    elif (mlat == lat-1):
        # 後邊界，故使用向後差分法
        vary= (ver[mlev,mlat,mlon] - ver[mlev,mlat-1,mlon]) / dy
        return vary
    else:
        # 其他非邊界項，使用中差法
        vary= (ver[mlev,mlat+1,mlon] - ver[mlev,mlat-1,mlon]) / (2*dy)
        return vary

#讀output.bin(二進位)
df=pd.DataFrame(np.fromfile('output.bin', np.dtype([('variable', 'float32')])))

#寫空的矩陣資料以便後面進行填寫
h = np.zeros((lev,lat,lon))
u = np.zeros((lev,lat,lon))
v = np.zeros((lev,lat,lon))
t = np.zeros((lev,lat,lon))

#寫資料入矩陣
a = 0
for mlev in range(lev):
    for mlat in range(lat):
        for mlon in range(lon):
            h[mlev,mlat,mlon] = df['variable'][a]
            u[mlev,mlat,mlon] = df['variable'][a+6125]    # 第6126個為第二個參數，6125=25*49*5
            v[mlev,mlat,mlon] = df['variable'][a+6125*2]
            t[mlev,mlat,mlon] = df['variable'][a+6125*3]
            a += 1

#寫空的散度矩陣資料以便後面進行填寫
div  = np.zeros((lev,lat,lon)) #散度

for mlev in range(lev):
    for mlat in range(lat):
        for mlon in range(lon):
            # 散度
            div[mlev,mlat,mlon]  = ddx(u,mlev,mlat,mlon) + ddy(v,mlev,mlat,mlon)

p=[1010,925,775,600,400,100] #垂直速度的氣壓面
#取120E的散度(16=120E)
lon120 = int((120 -90)/1.875)
div120=div[::,::,lon120]
#寫空的垂直速度矩陣資料以便後面進行填寫
w=np.zeros((len(p),lat))

for mlev in range(len(p)):
    for mlat in range(lat):
        if ( mlev ==0 ):
            w[mlev,mlat] = 0. #設定W0=0
            continue
        w[mlev,mlat] = w[mlev-1,mlat] + div120[mlev-1,mlat] * (p[mlev]-p[mlev-1])

#有初始風場後開始修正
e = np.zeros(lat)
dp = p[-1] - p[0]

# 修正 Divergence
for mlat in range(lat):
    # 修正參數
    e[mlat] = w[-1,mlat]/dp
    # 修正散度
    div[:,mlat] += -e[mlat]

# 修正 W
for mlev in range(len(p)):
    for mlat in range(lat):
        if ( mlev ==0 ):
            w[mlev,mlat] = 0. #設定W0=0
            continue
        w[mlev,mlat] = w[mlev-1,mlat] + div120[mlev-1,mlat] * (p[mlev]-p[mlev-1])

#繪圖
y = np.arange(lat)
z = p
Y,Z=np.meshgrid(y,z)
z1 = range(len(p))
Y,Z1 = np.meshgrid(y,z1)

plt.figure(figsize=(16,9))
plt.contour(Y*1.875+15,Z,w[Z1,Y],9,cmap = 'hsv') #hsv:彩虹color
plt.title('Vertical Velocity(120E)',fontsize=20)
plt.gca().invert_yaxis()
plt.yticks([1000,850,700,500,100])
plt.savefig('Vertical Velocity.jpg')
