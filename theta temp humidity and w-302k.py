import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 

filename = 'RCE_AT_302K.nc'

nc = netCDF4.Dataset(filename)
print nc.variables["dt_vertadv"]

senht = nc.variables["sens_ht"]
latht = nc.variables["lat_ht"]
tot_pr =  nc.variables["tot_precip"]

w = nc.variables["w"]
temp = nc.variables["T"]
qv = nc.variables["q"]
hei_th = nc.variables["h_rho"]
pre_th = nc.variables["p_theta"]
theta = nc.variables["theta"]

nt = 7200
nz = 63

pt_zt = np.zeros([63,7200])
sh_zt = np.zeros(7200)
lh_zt = np.zeros(7200)
t = np.zeros(7200)
temp_zt = np.zeros([63,7200])
qv_zt = np.zeros([63,7200])
hei_th_zt = np.zeros([63,7200])
pre_th_zt = np.zeros([63,7200])
w_zt = np.zeros([63,7200])
tot_pr_zt = np.zeros(7200)

for j in range(0,nt):
    sh_zt[j] = senht[j,0,0,0]
    lh_zt[j] = latht[j,0,0,0]
    t[j] = j
    for i in range(0,nz):
        pt_zt[i,j] = theta[j,i,0,0]
        temp_zt[i,j] = temp[j,i,0,0]
        qv_zt[i,j] = qv[j,i,0,0]
        hei_th_zt[i,j] = hei_th[j,i,0,0]
        pre_th_zt[i,j] = pre_th[j,i,0,0]
        w_zt[i,j] = w[j,i,0,0]
        
        
pre_th_mean = np.zeros(63)
for i in range(0,nz):
    pre_th_mean[i] = np.mean(pre_th_zt[i,:])
    
plt.figure(figsize=(7,10), dpi=80)
p1 = plt.subplot(221)
p1.plot(pt_zt[0:nz-5,4320],pre_th_mean[0:nz-5]/100.,linewidth=2,color='black')
p1.set_ylim(ymin=1000,ymax=0)
p1.set_xlim(xmin=290,xmax=450)
xmajorLocator=MultipleLocator(50)
p1.xaxis.set_major_locator(xmajorLocator)  
p1.set_xlabel('theta(K)')
p1.set_ylabel('Height(hPa)')
p2 = plt.subplot(222)
p2.plot(1000*qv_zt[0:nz-5,4320],pre_th_mean[0:nz-5]/100.,linewidth=2,color='black')
p2.set_ylim(ymin=1000,ymax=0)
p2.set_xlim(xmin=0,xmax=22)
p2.set_xlabel('Humidity(g/kg)')
p3 = plt.subplot(223)
p3.plot(temp_zt[0:nz-5,4320],pre_th_mean[0:nz-5]/100.,linewidth=2,color='black')
p3.set_ylim(ymin=1000,ymax=0)
p3.set_xlim(xmin=130,xmax=310)
xmajorLocator=MultipleLocator(50)
p3.xaxis.set_major_locator(xmajorLocator)  
p3.set_xlabel('Temperature(K)')
p3.set_ylabel('Height(hPa)')
p4 = plt.subplot(224)
p4.plot(w_zt[0:nz-5,4320],pre_th_mean[0:nz-5]/100.,linewidth=2,color='black')
p4.set_ylim(ymin=1000,ymax=0)
p4.set_xlabel('Vertical velocity(m/s)')
plt.savefig('theta temp humidity and w-302k.pdf')
plt.show()




