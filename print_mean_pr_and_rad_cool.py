import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 

filename = '../298_R298.nc'

nc = netCDF4.Dataset(filename)
rain = nc.variables["ls_rain3d"]
print rain
nt = 7200
nz = 63

#rain_zt = np.zeros([nz,nt])
#for j in range(0,nt):
#    for i in range(0,nz):
#        rain_zt[i,j] = rain[j,i,0,0]
#        
#        
#        
#        
#senht = nc.variables["sens_ht"]
#latht = nc.variables["lat_ht"]
tot_pr =  nc.variables["tot_precip"]
#
#
#sh_zt = np.zeros(nt)
#lh_zt = np.zeros(nt)
#t = np.zeros(nt)
tot_pr_zt = np.zeros(nt)
#
for j in range(0,nt):
    tot_pr_zt[j] = tot_pr[j,0,0,0]
#    sh_zt[j] = senht[j,0,0,0]
#    lh_zt[j] = latht[j,0,0,0]
#    t[j] = j
#
#plt.figure(figsize=(7,4.5),dpi=80)
#plt.plot(t/144.,tot_pr_zt[:],linewidth=1,color='blue')
#plt.xlim([0,50])
#plt.xlabel('Time(days)')
#plt.ylabel('Total Precipitation Rate(kg/m2/day)')
#plt.savefig('total precipitation rate(50d simulation).pdf')
#plt.show()
#
#
#plt.figure(figsize=(7,4.5),dpi=80)
#plt.plot(t/144.,sh_zt[:],linewidth=1,color='blue',label='sensible heaing')
#plt.plot(t/144.,lh_zt[:],linewidth=1,color='red',label='latent heaing')
#plt.xlim([0,50])
#plt.xlabel('Time(days)')
#plt.ylabel('Sensible and Laten Heating(W/m2)')
#plt.legend()
#plt.savefig('Suface heating(50d simulation).pdf')
#plt.show()
#
totforc = nc.variables["dt_totforc"]
pre_th = nc.variables["p_theta"]
sqwtg = nc.variables["dq_vertadv"]

neq = 4320
dt_scm = 600.
gra = 9.8
L_vap = 2.501e6    
C_p = 1005.0     
 
t = np.zeros(7200)
totforc_zt = np.zeros([nz,nt])
sqwtg_zt = np.zeros([nz,nt])
totforc_mean = np.zeros(nz)
pre_th_zt = np.zeros([nz,nt])
pre_th_mean = np.zeros(nz)
sqwtg_mean = np.zeros(nz)
dp_mean = np.zeros(nz)
dp = np.zeros([nz,nt])

ppt_tot = np.zeros(nt)
#
for j in range(0,nt):
    t[j] = j
    for i in range(0,nz):
        totforc_zt[i,j] = totforc[j,i,0,0]
        sqwtg_zt[i,j] = sqwtg[j,i,0,0]
        pre_th_zt[i,j] = pre_th[j,i,0,0]
       
for i in range(0,nz):
    totforc_mean[i] = np.mean(totforc_zt[i,neq:nt])/dt_scm
    sqwtg_mean[i] = np.mean(sqwtg_zt[i,neq:nt])/dt_scm
    pre_th_mean[i] = np.mean(pre_th_zt[i,neq:nt])

for i in range(0,nz-1):        
    dp_mean[i+1] = -(pre_th_mean[i+1] - pre_th_mean[i])
 
rad_cool = 0.
moist_wtg = 0.   
for i in range(0,nz-1):
    rad_cool = rad_cool + (totforc_mean[i+1] + totforc_mean[i])*dp_mean[i+1]/2.
    moist_wtg = moist_wtg + (sqwtg_mean[i+1] + sqwtg_mean[i])*dp_mean[i+1]/2.
    
rad_cool = -86400.*(C_p/L_vap)*rad_cool/gra
moist_wtg = 86400*moist_wtg/gra
mean_pr = np.mean(tot_pr_zt[neq:nt])

for j in range(0,nt):    
#    ppt_tot[j] = 0.
    for i in range(0,nz-1):
        dp[i+1,j] = -(pre_th_zt[i+1,j] - pre_th_zt[i,j])
#        ppt_tot[j] = ppt_tot[j] + (dqconvertp_zt[i,j] + dqconvertp_zt[i+1,j])*dp[i+1,j]/2.
#
#ppt_mean = -86400*np.mean(ppt_tot[neq:nt])/gra
#
print ('The radiative cooling is', rad_cool)
print ('The mean precipitation is', mean_pr)
#print ('The mean converted precipitaiton is', ppt_mean)    

#rc = np.zeros(nt)
#for j in range(0,nt):
#    rc[j] = 0.
#    for i in range(0,nz-1):
#        rc[j] = rc[j] + (totforc_zt[i+1,j] + totforc_zt[i,j])*dp[i+1,j]/2.
#    rc[j] = -144*(C_p/L_vap)*rc[j]/gra
#
#tot_pr_3h = np.zeros(50)
#time = np.zeros(50)
#for i in range(0,50):
#    time[i] = i
#    a = i*144
#    tot_pr_3h[i] = np.mean(tot_pr_zt[a:a+143])
#    
#plt.figure(figsize=(7,4.5),dpi=80)
#plt.plot(time,tot_pr_3h[:],linewidth=1,color='blue',label='precipitation rate at 298K')
#plt.plot(t/144.,rc[:],linewidth=2,color='red',label='radiative cooling')
#plt.xlim([0,50])
#plt.ylim([0,7])
#plt.legend(bbox_to_anchor=(1.0, 0.2))
#plt.xlabel('Time(days)')
#plt.ylabel('Rain Rates(mm/day)')
#plt.savefig('total precipitation rate(50d simulation)-1d mean(298K).pdf')
#plt.show()
#
#
