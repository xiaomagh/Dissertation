import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 

filename = 'gabls2_murkem_tracer_l63.nc'

nc = netCDF4.Dataset(filename)

temp = nc.variables["T"]
qv = nc.variables["q"]
hei_th = nc.variables["h_rho"]
pre_th = nc.variables["p_theta"]
theta = nc.variables["theta"]
rh = nc.variables["rh2"]

nt = 7200
nz = 63
neq = 4320

pt_zt = np.zeros([nz,nt])
t = np.zeros(nt)
temp_zt = np.zeros([nz,nt])
qv_zt = np.zeros([nz,nt])
hei_th_zt = np.zeros([nz,nt])
pre_th_zt = np.zeros([nz,nt])
rh_zt = np.zeros([nz,nt])

for j in range(0,nt):
    t[j] = j
    for i in range(0,nz):
        pt_zt[i,j] = theta[j,i,0,0]
        temp_zt[i,j] = temp[j,i,0,0]
        qv_zt[i,j] = qv[j,i,0,0]
        hei_th_zt[i,j] = hei_th[j,i,0,0]
        pre_th_zt[i,j] = pre_th[j,i,0,0]   
        rh_zt[i,j] = rh[j,i,0,0]   
        
pre_th_mean = np.zeros(63)
pt_mean = np.zeros(63)
temp_mean = np.zeros(63)
q_mean = np.zeros(63)
rh_mean = np.zeros(63)
for i in range(0,nz):
    pre_th_mean[i] = np.mean(pre_th_zt[i,neq:nt])
    pt_mean[i] = np.mean(pt_zt[i,neq:nt])
    temp_mean[i] = np.mean(temp_zt[i,neq:nt])
    q_mean[i] = np.mean(qv_zt[i,neq:nt])
    rh_mean[i] = np.mean(rh_zt[i,neq:nt])

#def e_saturated_w(t):
#    es = 6.112*np.exp(17.67*(t-273.15)/(t-273.15+243.5))*100
#    return es
#
#def e_saturated_i(t):
#    es = 6.112*10**(9.5*(t-273.15)/(t-273.15+265.5))*100
#    return es
#
#for i in range(0,nz):
#    es = e_saturated_w(temp_mean[i])
#    ws = 0.62197*(es/(pre_th_mean[i]-es))
#    rh_mean[i] = 100*q_mean[i]/ws
#for i in range(nz,nz):
#    es = e_saturated_i(temp_mean[i])
#    ws = 0.62197*(es/(pre_th_mean[i]-es))
#    rh_mean[i] = 100*q_mean[i]/ws
#   
plt.figure(figsize=(7,10), dpi=80)
p1 = plt.subplot(221)
p1.plot(pt_mean[0:nz-5],pre_th_mean[0:nz-5]/100.,linewidth=2,color='black')
p1.set_ylim(ymin=1000,ymax=0)
p1.set_xlim(xmin=290,xmax=450)
xmajorLocator=MultipleLocator(50)
p1.xaxis.set_major_locator(xmajorLocator)  
p1.set_xlabel('theta(K)')
p1.set_ylabel('Height(hPa)')
p2 = plt.subplot(222)
p2.plot(1000*q_mean[0:nz-5],pre_th_mean[0:nz-5]/100.,linewidth=2,color='black')
p2.set_ylim(ymin=1000,ymax=0)
p2.set_xlim(xmin=0,xmax=22)
p2.set_xlabel('Humidity(g/kg)')
p3 = plt.subplot(223)
p3.plot(temp_mean[0:nz-5],pre_th_mean[0:nz-5]/100.,linewidth=2,color='black')
p3.set_ylim(ymin=1000,ymax=0)
p3.set_xlim(xmin=130,xmax=310)
xmajorLocator=MultipleLocator(50)
p3.xaxis.set_major_locator(xmajorLocator)  
p3.set_xlabel('Temperature(K)')
p3.set_ylabel('Height(hPa)')
p4 = plt.subplot(224)
p4.plot(100*rh_mean[0:nz-5],pre_th_mean[0:nz-5]/100.,linewidth=2,color='black')
p4.set_ylim(ymin=1000,ymax=0)
p4.set_xlim(xmin=0,xmax=100)
p4.set_xlabel('Relative Humidity(%)')
plt.savefig('the average theta temp and humidity of last 20d.pdf')
plt.show()




