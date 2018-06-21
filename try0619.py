import netCDF4
import numpy as np
import matplotlib.pyplot as plt

filename = 'gabls2_murkem_tracer_l63.nc'

nc = netCDF4.Dataset(filename)
print nc.variables["rh2"]
c = nc.variables["rh2"]
clo = np.zeros([63,7200])

for j in range(0,7200):
    
    for i in range(0,63):
        clo[i,j] = c[j,i,0,0]
#dcflpc2 = nc.variables["dcfl_pc2ck"]
#dcffpc2 = nc.variables["dcff_pc2ck"]
#dcflconv = nc.variables["dcfl_conv"]
#dcffconv = nc.variables["dcff_conv"]
#dcfllsr = nc.variables["dcfl_lsr"]
#dcfflsr = nc.variables["dcff_lsr"]
#pre_th = nc.variables["p_theta"]
#
#nt = 7200
#nz = 63
#neq = 4320
#
#t = np.zeros(nt)
#dcfl = np.zeros([nz,nt])
#dcff = np.zeros([nz,nt])
#pre_th_zt = np.zeros([nz,nt])
#dcflpc2_zt = np.zeros([nz,nt])
#dcffpc2_zt = np.zeros([nz,nt])
#dcflconv_zt = np.zeros([nz,nt])
#dcffconv_zt = np.zeros([nz,nt])
#dcfllsr_zt = np.zeros([nz,nt])
#dcfflsr_zt = np.zeros([nz,nt])
#
#
#for j in range(0,nt):
#    t[j] = j
#    for i in range(0,nz):
#        pre_th_zt[i,j] = pre_th[j,i,0,0]
#        dcflpc2_zt[i,j] = dcflpc2[j,i,0,0]
#        dcffpc2_zt[i,j] = dcffpc2[j,i,0,0]
#        dcflconv_zt[i,j] = dcflconv[j,i,0,0]
#        dcffconv_zt[i,j] = dcffconv[j,i,0,0]
#        dcfllsr_zt[i,j] = dcfllsr[j,i,0,0]
#        dcfflsr_zt[i,j] = dcfflsr[j,i,0,0]
#        dcfl[i,j] = dcflpc2_zt[i,j] + dcflconv_zt[i,j] + dcfllsr_zt[i,j]
#        dcff[i,j] = dcffpc2_zt[i,j] + dcffconv_zt[i,j] + dcfflsr_zt[i,j]
#
#pre_th_mean = np.zeros(63)
#for i in range(0,nz):
#    pre_th_mean[i] = np.mean(pre_th_zt[i,neq:nt])
#
#plt.figure(figsize=(8,7), dpi=80)
#p1 = plt.subplot(211)
#et = p1.contourf(t[:]/144.,pre_th_mean[0:nz-5]/100.,dcfl[0:nz-5,:],15,
#                 vmax = 0.75, vmin = 0, cmap='OrRd')
#plt.colorbar(et,extend='both')
#p1.set_ylim(ymin=1000,ymax=0)
#p1.set_xlabel('Time(day)')
#p1.set_ylabel('Height(hPa)')
#p2 = plt.subplot(212)
#ft = p2.contourf(t[:]/144.,pre_th_mean[0:nz-5]/100.,dcff[0:nz-5,:],15,
#                 vmax = 0.6, vmin = 0, cmap='OrRd')
#plt.colorbar(ft,extend='both')
#p2.set_ylim(ymin=1000,ymax=0)
#p2.set_xlabel('Time(day)')
#p2.set_ylabel('Height(hPa)')
#plt.savefig('contour of liquid and frozen cloud.pdf')
#plt.show()




