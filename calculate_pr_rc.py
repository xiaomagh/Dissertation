import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 

filename = 'gabls2_murkem_tracer_l63-5d.nc'

nc = netCDF4.Dataset(filename)

totforc = nc.variables["dt_totforc"]
pre_th = nc.variables["p_theta"]
tot_pr =  nc.variables["tot_precip"]
dqclpc2 = nc.variables["dqcl_pc2ck"]
dqcfpc2 = nc.variables["dqcf_pc2ck"]
dqpc2 = nc.variables["dq_pc2ck"]
dqclconv = nc.variables["dqcl_conv"]
dqcfconv = nc.variables["dqcf_conv"]
dqconv = nc.variables["dq_conv"]
dqcllsr = nc.variables["dqcl_lsr"]
dqcflsr = nc.variables["dqcf_lsr"]
dqlsr = nc.variables["dq_lsr"]

nt = 720
nz = 63
neq = 620
dt_scm = 600.
gra = 9.8
L_vap = 2.501e6    
C_p = 1005.0      

totforc_zt = np.zeros([63,720])
t = np.zeros(720)
totforc_mean = np.zeros(63)
pre_th_zt = np.zeros([63,720])
pre_th_mean = np.zeros(63)
dp_mean = np.zeros(63)
dp = np.zeros([63,720])
tot_pr_zt = np.zeros(720)
dqclpc2_zt = np.zeros([63,720])
dqcfpc2_zt = np.zeros([63,720])
dqpc2_zt = np.zeros([63,720])
dqclconv_zt = np.zeros([63,720])
dqcfconv_zt = np.zeros([63,720])
dqconv_zt = np.zeros([63,720])
dqcllsr_zt = np.zeros([63,720])
dqcflsr_zt = np.zeros([63,720])
dqlsr_zt = np.zeros([63,720])
ppt_tot = np.zeros(720)

for j in range(0,nt):
    tot_pr_zt[j] = tot_pr[j,0,0,0]
    t[j] = j
#    for i in range(0,nz):
#        totforc_zt[i,j] = totforc[j,i,0,0]
#        pre_th_zt[i,j] = pre_th[j,i,0,0]
#        dqclpc2_zt[i,j] = dqclpc2[j,i,0,0]
#        dqcfpc2_zt[i,j] = dqcfpc2[j,i,0,0]
#        dqpc2_zt[i,j] = dqpc2[j,i,0,0]
#        dqclconv_zt[i,j] = dqclconv[j,i,0,0]
#        dqcfconv_zt[i,j] = dqcfconv[j,i,0,0]
#        dqconv_zt[i,j] = dqconv[j,i,0,0]
#        dqcllsr_zt[i,j] = dqcllsr[j,i,0,0]
#        dqcflsr_zt[i,j] = dqcflsr[j,i,0,0]
#        dqlsr_zt[i,j] = dqlsr[j,i,0,0]
#        
#for i in range(0,nz):
#    totforc_mean[i] = np.mean(totforc_zt[i,neq:nt])/dt_scm
#    pre_th_mean[i] = np.mean(pre_th_zt[i,neq:nt])
#
#for i in range(0,nz-1):        
#    dp_mean[i+1] = -(pre_th_mean[i+1] - pre_th_mean[i])
# 
#rad_cool = 0.   
#for i in range(0,nz-1):
#    rad_cool = rad_cool + (totforc_mean[i+1] + totforc_mean[i])*dp_mean[i+1]/2.
#    
#rad_cool = -86400.*(C_p/L_vap)*rad_cool/gra
#mean_pr = np.mean(tot_pr_zt[neq:nt])
#
#dqconvertp_zt = np.zeros([63,720])
#for j in range(0,nt):
#    for i in range(0,nz):
#        dqconvertp_zt[i,j] = (dqclpc2_zt[i,j]+dqcfpc2_zt[i,j]+dqpc2_zt[i,j]+  \
#                            dqclconv_zt[i,j]+dqcfconv_zt[i,j]+dqconv_zt[i,j]+ \
#                            dqcllsr_zt[i,j]+dqcflsr_zt[i,j]+dqlsr_zt[i,j])/dt_scm
#where_are_nan = np.isnan(dqconvertp_zt)
#dqconvertp_zt[where_are_nan] = 0.  
#
#for j in range(0,nt):    
#    ppt_tot[j] = 0.
#    for i in range(0,nz-1):
#        dp[i+1,j] = -(pre_th_zt[i+1,j] - pre_th_zt[i,j])
#        ppt_tot[j] = ppt_tot[j] + (dqconvertp_zt[i,j] + dqconvertp_zt[i+1,j])*dp[i+1,j]/2.
#
#ppt_mean = -86400*np.mean(ppt_tot[neq:nt])/gra
#
#print ('The radiative cooling is', rad_cool)
#print ('The mean precipitation is', mean_pr)
#print ('The mean converted precipitaiton is', ppt_mean)    
#
#rc = np.zeros(720)
#for j in range(0,nt):
#    rc[j] = 0.
#    for i in range(0,nz-1):
#        rc[j] = rc[j] + (totforc_zt[i+1,j] + totforc_zt[i,j])*dp[i+1,j]/2.
#    rc[j] = -144*(C_p/L_vap)*rc[j]/gra

tot_pr_3h = np.zeros(40)
time = np.zeros(40)
for i in range(0,40):
    time[i] = i
    a = i*18
    tot_pr_3h[i] = np.mean(tot_pr_zt[a:a+17])
    
plt.figure(figsize=(7,4.5),dpi=80)
plt.plot(time/8.,tot_pr_3h[:],linewidth=2,color='blue',label='precipitation')
plt.xlim([0,5])
plt.xlabel('Time(days)')
plt.ylabel('Total Precipitation Rate(kg/m2/day)')
plt.savefig('total precipitation rate.pdf')
plt.show()



    