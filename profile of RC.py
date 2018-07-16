import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 

filename1 = 'RCE_AT_298K.nc'
filename2 = 'RCE_AT_300K.nc'
filename3 = 'RCE_AT_302K.nc'

nc1 = netCDF4.Dataset(filename1)
nc2 = netCDF4.Dataset(filename2)
nc3 = netCDF4.Dataset(filename3)

rc1 = nc1.variables["dt_totforc"]
rc2 = nc2.variables["dt_totforc"]
rc3 = nc3.variables["dt_totforc"]
pre_th1 = nc1.variables["p_theta"]
pre_th2 = nc2.variables["p_theta"]
pre_th3 = nc3.variables["p_theta"]
nt = 7200
nz = 63
rc1_zt = np.zeros([nz,nt])
rc2_zt = np.zeros([nz,nt])
rc3_zt = np.zeros([nz,nt])
pre_th_zt1 = np.zeros([63,7200])
pre_th_zt2 = np.zeros([63,7200])
pre_th_zt3 = np.zeros([63,7200])
for j in range(0,nt):
    for i in range(0,nz):
        rc1_zt[i,j] = rc1[j,i,0,0]
        rc2_zt[i,j] = rc2[j,i,0,0]
        rc3_zt[i,j] = rc3[j,i,0,0]
        pre_th_zt1[i,j] = pre_th1[j,i,0,0]
        pre_th_zt2[i,j] = pre_th2[j,i,0,0]
        pre_th_zt3[i,j] = pre_th3[j,i,0,0]
        
pre_th_mean1 = np.zeros(63)
pre_th_mean2 = np.zeros(63)
pre_th_mean3 = np.zeros(63)
rc1_mean = np.zeros(63)
rc2_mean = np.zeros(63)
rc3_mean = np.zeros(63)
for i in range(0,nz):
    pre_th_mean1[i] = np.mean(pre_th_zt1[i,:])
    pre_th_mean2[i] = np.mean(pre_th_zt2[i,:])
    pre_th_mean3[i] = np.mean(pre_th_zt3[i,:])
    rc1_mean[i] = np.mean(rc1_zt[i,4320:nt])
    rc2_mean[i] = np.mean(rc2_zt[i,4320:nt])    
    rc3_mean[i] = np.mean(rc3_zt[i,4320:nt])
    
for i in range(0,nz):
    rc1_mean[i] = 86400*rc1_mean[i]/600. 
    rc2_mean[i] = 86400*rc2_mean[i]/600. 
    rc3_mean[i] = 86400*rc3_mean[i]/600. 
    
plt.figure(figsize=(3.5,5), dpi=80)
plt.plot(rc1_mean[0:nz-5],pre_th_mean1[0:nz-5]/100.,linewidth=2,color='black',label='298K')
plt.plot(rc2_mean[0:nz-5],pre_th_mean2[0:nz-5]/100.,linewidth=2,color='blue',label='300K')
plt.plot(rc3_mean[0:nz-5],pre_th_mean3[0:nz-5]/100.,linewidth=2,color='red',label='302K')
plt.axvline(0, linestyle=':', color='black')
plt.ylim([1000,0])
plt.xlim([-1.6,1.6])  
plt.legend(bbox_to_anchor=(1.0, 0.4))
plt.xlabel('Radiative Cooling(K/day)')
plt.ylabel('Height(hPa)')
plt.savefig('profile of RC.pdf')
plt.show()

    
    