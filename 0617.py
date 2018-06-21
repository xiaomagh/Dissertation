import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 

filename = 'gabls2_murkem_tracer_l63.nc'

nc = netCDF4.Dataset(filename)

w = nc.variables["w"]
temp = nc.variables["T"]
qv = nc.variables["q"]
hei_th = nc.variables["h_rho"]
pre_th = nc.variables["p_theta"]
theta = nc.variables["theta"]

nt = 7200
nz = 63
neq = 4320

pt_zt = np.zeros([nz,nt])
t = np.zeros(nt)
temp_zt = np.zeros([nz,nt])
hei_th_zt = np.zeros([nz,nt])
pre_th_zt = np.zeros([nz,nt])
theta_diff_zt = np.zeros([nz,nt])
temp_diff_zt = np.zeros([nz,nt])

for j in range(0,nt):
    t[j] = j
    for i in range(0,nz):
        pt_zt[i,j] = theta[j,i,0,0]
        pre_th_zt[i,j] = pre_th[j,i,0,0]
        temp_zt[i,j] = temp[j,i,0,0] 
        
where_are_nan = np.isnan(pre_th_zt)
pre_th_zt[where_are_nan] = 0.                
        
pre_th_mean = np.zeros(63)
pt_mean = np.zeros(63)
temp_mean = np.zeros(63)

for i in range(0,nz):
    pre_th_mean[i] = np.mean(pre_th_zt[i,neq:nt])
    pt_mean[i] = np.mean(pt_zt[i,neq:nt])
    temp_mean[i] = np.mean(temp_zt[i,neq:nt])
    for j in range(0,nt):
        theta_diff_zt[i,j] = pt_zt[i,j] - pt_mean[i]
        temp_diff_zt[i,j] = temp_zt[i,j] - temp_mean[i]


plt.figure(figsize=(8,7), dpi=80)
p1 = plt.subplot(211)
et = p1.contourf(t[:]/144.,pre_th_mean[:]/100.,theta_diff_zt, 20,
            vmin=-10, vmax=10, cmap='seismic')
plt.colorbar(et,extend='both')
p1.set_ylim(ymin=1000,ymax=0)
p1.set_xlabel('Time(day)')
p1.set_ylabel('Height(hPa)')
p2 = plt.subplot(212)
ft = p2.contourf(t[:]/144.,pre_th_mean[:]/100.,temp_diff_zt, 20,
            vmin=-10, vmax=10, cmap='seismic')
plt.colorbar(ft,extend='both')
p2.set_ylim(ymin=1000,ymax=0)
p2.set_xlabel('Time(day)')
p2.set_ylabel('Height(hPa)')
plt.savefig('contour of theta and temp diff.pdf')
plt.show()
