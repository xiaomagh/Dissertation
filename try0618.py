import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 

filename = 'gabls2_murkem_tracer_l63.nc'

nc = netCDF4.Dataset(filename)

qv = nc.variables["q"]
hei_th = nc.variables["h_rho"]
pre_th = nc.variables["p_theta"]
temp = nc.variables["T"]

nt = 7200
nz = 63
neq = 4320

t = np.zeros(nt)
temp_zt = np.zeros([nz,nt])
qv_zt = np.zeros([nz,nt])
hei_th_zt = np.zeros([nz,nt])
pre_th_zt = np.zeros([nz,nt])
rh_zt = np.zeros([nz,nt])
q_diff_zt = np.zeros([nz,nt])
rh_diff_zt = np.zeros([nz,nt])

for j in range(0,nt):
    t[j] = j
    for i in range(0,nz-5):
        temp_zt[i,j] = temp[j,i,0,0]
        qv_zt[i,j] = qv[j,i,0,0]
        hei_th_zt[i,j] = hei_th[j,i,0,0]
        pre_th_zt[i,j] = pre_th[j,i,0,0]     

def e_saturated_w(t):
    es = 6.112*np.exp(17.67*(t-273.15)/(t-273.15+243.5))*100
    return es

def e_saturated_i(t):
    es = 6.112*10**(9.5*(t-273.15)/(t-273.15+265.5))*100
    return es

for j in range(0,nt):
    for i in range(0,41):
        es = e_saturated_w(temp_zt[i,j])
        ws = 0.62197*(es/(pre_th_zt[i,j]-es))
        rh_zt[i,j] = 100*qv_zt[i,j]/ws
    for i in range(41,nz):
        es = e_saturated_i(temp_zt[i,j])
        ws = 0.62197*(es/(pre_th_zt[i,j]-es))
        rh_zt[i,j] = 100*qv_zt[i,j]/ws

    
pre_th_mean = np.zeros(63)
q_mean = np.zeros(63)
rh_mean = np.zeros(63)
for i in range(0,nz):
    pre_th_mean[i] = np.mean(pre_th_zt[i,neq:nt])
    q_mean[i] = np.mean(qv_zt[i,neq:nt])
    rh_mean[i] = np.mean(rh_zt[i,neq:nt])
    for j in range(0,nt):
        q_diff_zt[i,j] = qv_zt[i,j] - q_mean[i]
        rh_diff_zt[i,j] = rh_zt[i,j] - rh_mean[i]

plt.figure(figsize=(8,7), dpi=80)
p1 = plt.subplot(211)
et = p1.contourf(t[:]/144.,pre_th_mean[0:nz-5]/100.,q_diff_zt[0:nz-5,:]*1000, 20,
            vmin=-5, vmax=5, cmap='RdYlBu')
plt.colorbar(et,extend='both')
p1.set_ylim(ymin=1000,ymax=0)
p1.set_xlabel('Time(day)')
p1.set_ylabel('Height(hPa)')
p2 = plt.subplot(212)
ft = p2.contourf(t[:]/144.,pre_th_mean[0:nz-20]/100.,rh_diff_zt[0:nz-20,:], 20,
                 vmin=-50, vmax=50,cmap='RdYlBu')
plt.colorbar(ft,extend='both')
p2.set_ylim(ymin=1000,ymax=0)
p2.set_xlabel('Time(day)')
p2.set_ylabel('Height(hPa)')
plt.savefig('contour of q and rh diff.pdf')
plt.show()








