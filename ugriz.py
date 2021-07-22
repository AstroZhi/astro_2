import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import scipy.interpolate as spi
import math

curves = fits.open(r'C:\Users\DELL\Desktop\astro\filter_curves.fits')
star = fits.open(r'C:\Users\DELL\Desktop\astro\spec-55859-F5902_sp01-014.fits')

#提取元素
u = curves[1].data
g = curves[2].data
r = curves[3].data
i = curves[4].data
z = curves[5].data

flux_g = g.field(1)
wavelength_g = g.field(0)

flux_r = r.field(1)
wavelength_r = r.field(0)

wavelength_star = np.zeros((1, 3908))
flux_star = np.zeros((1, 3908))

wavelength_star = star[0].data[2]
print(wavelength_star)
flux_star = star[0].data[0]
print(flux_star)

#插值
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

# 数据准备
print(wavelength_g.max())
print(wavelength_g.min())
print(wavelength_r.max())
print(wavelength_r.min())

print(wavelength_star.max())
print(wavelength_star.min())

new_xg = np.arange(3630, 5830, 0.01)
new_xr = np.arange(5380, 7230, 0.01)
new_xstar = np.arange(3699, 9100, 0.01)

#进行一阶样条插值
ipo1 = spi.splrep(wavelength_star, flux_star, k=1) #样本点导入，生成参数
new_ystar = spi.splev(new_xstar, ipo1) #根据观测点和样条参数，生成插值

#进行三次样条插值
ipo_g = spi.splrep(wavelength_g, flux_g, k=3) #样本点导入，生成参数
iy_g = spi.splev(new_xg, ipo_g) #根据观测点和样条参数，生成插值
ipo_r = spi.splrep(wavelength_r, flux_r, k=3) #样本点导入，生成参数
iy_r = spi.splev(new_xr, ipo_r) #根据观测点和样条参数，生成插值

g = 0
r = 0
for i in range(213100):
    ind = i+6900
    g += new_ystar[i] * iy_g[ind] * 0.01
print(g)
for i in range(185000):
    ind = i + 168100
    r += new_ystar[ind] * iy_r[i] * 0.01
print(r)

mag = 2.5*math.log10(r/g)
print(mag)






