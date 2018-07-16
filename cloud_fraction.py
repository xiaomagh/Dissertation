import netCDF4
import numpy as np
import matplotlib.pyplot as plt

filename = 'gabls2_murkem_tracer_l63.nc'

nc = netCDF4.Dataset(filename)

dcflpc2 = nc.variables["dcfl_pc2ck"]
dcffpc2 = nc.variables["dcff_pc2ck"]
dcflconv = nc.variables["dcfl_conv"]
dcffconv = nc.variables["dcff_conv"]
dcfllsr = nc.variables["dcfl_lsr"]
dcfflsr = nc.variables["dcff_lsr"]
pre_th = nc.variables["p_theta"]

nt = 7200
nz = 63
neq = 4320

t = np.zeros(nt)
dcfl = np.zeros([nz,nt])
dcff = np.zeros([nz,nt])
pre_th_zt = np.zeros([nz,nt])
dcflpc2_zt = np.zeros([nz,nt])
dcffpc2_zt = np.zeros([nz,nt])
dcflconv_zt = np.zeros([nz,nt])
dcffconv_zt = np.zeros([nz,nt])
dcfllsr_zt = np.zeros([nz,nt])
dcfflsr_zt = np.zeros([nz,nt])


for j in range(0,nt):
    t[j] = j
    for i in range(0,nz):
        pre_th_zt[i,j] = pre_th[j,i,0,0]
        dcflpc2_zt[i,j] = dcflpc2[j,i,0,0]
        dcffpc2_zt[i,j] = dcffpc2[j,i,0,0]
        dcflconv_zt[i,j] = dcflconv[j,i,0,0]
        dcffconv_zt[i,j] = dcffconv[j,i,0,0]
        dcfllsr_zt[i,j] = dcfllsr[j,i,0,0]
        dcfflsr_zt[i,j] = dcfflsr[j,i,0,0]
        dcfl[i,j] = dcflpc2_zt[i,j] + dcflconv_zt[i,j] + dcfllsr_zt[i,j]
        if dcfl[i,j] < 0.:
            dcfl[i,j] = 0.
        dcff[i,j] = dcffpc2_zt[i,j] + dcffconv_zt[i,j] + dcfflsr_zt[i,j]
        if dcff[i,j] < 0.:
            dcff[i,j] = 0.

pre_th_mean = np.zeros(63)
dcfl_mean = np.zeros(63)
dcff_mean = np.zeros(63)
for i in range(0,nz):
    pre_th_mean[i] = np.mean(pre_th_zt[i,neq:nt])
    dcfl_mean[i] = np.mean(dcfl[i,neq:nt])
    dcff_mean[i] = np.mean(dcff[i,neq:nt])

plt.figure()
plt.plot(dcfl_mean[0:nz-5],pre_th_mean[0:nz-5]/100.,linewidth=2,color='blue',label='liquid cloud fraction')
plt.plot(dcff_mean[0:nz-5],pre_th_mean[0:nz-5]/100.,linewidth=2,color='red',label='frozen cloud fraction')
plt.ylim([1000,0])
plt.legend()
plt.savefig('profile of liquid and frozen cloud.pdf')
plt.show()




