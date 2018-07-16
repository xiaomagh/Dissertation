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

qv1 = nc1.variables["theta"]
qv2 = nc2.variables["theta"]
qv3 = nc3.variables["theta"]
pre_th1 = nc1.variables["p_theta"]
pre_th2 = nc2.variables["p_theta"]
pre_th3 = nc3.variables["p_theta"]
nt = 7200
nz = 63
q1_zt = np.zeros([nz,nt])
q2_zt = np.zeros([nz,nt])
q3_zt = np.zeros([nz,nt])
pre_th_zt1 = np.zeros([63,7200])
pre_th_zt2 = np.zeros([63,7200])
pre_th_zt3 = np.zeros([63,7200])
for j in range(0,nt):
    for i in range(0,nz):
        q1_zt[i,j] = qv1[j,i,0,0]
        q2_zt[i,j] = qv2[j,i,0,0]
        q3_zt[i,j] = qv3[j,i,0,0]
        pre_th_zt1[i,j] = pre_th1[j,i,0,0]
        pre_th_zt2[i,j] = pre_th2[j,i,0,0]
        pre_th_zt3[i,j] = pre_th3[j,i,0,0]
        
pre_th_mean1 = np.zeros(63)
pre_th_mean2 = np.zeros(63)
pre_th_mean3 = np.zeros(63)
q1_mean = np.zeros(63)
q2_mean = np.zeros(63)
q3_mean = np.zeros(63)
for i in range(0,nz):
    pre_th_mean1[i] = np.mean(pre_th_zt1[i,:])
    pre_th_mean2[i] = np.mean(pre_th_zt2[i,:])
    pre_th_mean3[i] = np.mean(pre_th_zt3[i,:])
    q1_mean[i] = np.mean(q1_zt[i,4320:nt])
    q2_mean[i] = np.mean(q2_zt[i,4320:nt])    
    q3_mean[i] = np.mean(q3_zt[i,4320:nt])
    
plt.figure(figsize=(3.5,5), dpi=80)
plt.plot(q1_mean[0:nz-5],pre_th_mean1[0:nz-5]/100.,linewidth=2,color='black',label='298K')
plt.plot(q2_mean[0:nz-5],pre_th_mean2[0:nz-5]/100.,linewidth=2,color='blue',label='300K')
plt.plot(q3_mean[0:nz-5],pre_th_mean3[0:nz-5]/100.,linewidth=2,color='red',label='302K')
plt.xlim([290,450])
plt.ylim([1000,0])
plt.legend(bbox_to_anchor=(1.0, 0.4))
plt.xlabel('Potential Temperature(K)')
plt.ylabel('Height(hPa)')
plt.savefig('profile of theta with diff SST.pdf')
plt.show()



