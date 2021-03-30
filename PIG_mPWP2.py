import numpy as np, matplotlib.pyplot as plt, seawater as sw, mPWP as PWP, time;
from datetime import datetime

# set up initial parameters (these are also set as defaults in the PWP code, so don't need to put them in)
T0 = 0; # reference temperature (degC)
S0 = 34; # reference salinity (parts per thousand; ppt)
rho0 = 1025; # reference density (kg m^-3)
alpha = -sw.alpha(S0,T0,0)*rho0 # thermal expansion   coefficient (kg m^-3 degC^-1)
beta = sw.beta(S0,T0,0)*rho0 # haline contraction coefficient (kg m^-3 ppt^-1 )
latitude = -75; # degrees north
f = sw.f(latitude) # planetary vorticity

fwflux1 = np.loadtxt('/home/erobo/ThompsonResearch/Data/fwflux1_oneyear.csv',delimiter=',')
qnet1 = np.loadtxt('/home/erobo/ThompsonResearch/Data/qnet1_oneyear.csv',delimiter=',')
taux1 = np.loadtxt('/home/erobo/ThompsonResearch/Data/taux1_oneyear.csv',delimiter=',')
tauy1 = np.loadtxt('/home/erobo/ThompsonResearch/Data/tauy1_oneyear.csv',delimiter=',')
temp_profile = np.loadtxt('/home/erobo/ThompsonResearch/Data/temp_profile_snapshot.txt');
salt_profile = np.loadtxt('/home/erobo/ThompsonResearch/Data/salt_profile_snapshot.txt');
depth_profile = np.loadtxt('/home/erobo/ThompsonResearch/Data/Schodlok_depths.txt');

dt = 1; # seconds
days = 365;
t = np.arange(0,86400*days,dt)

fwflux_interp =  np.interp(t/86400, range(364), fwflux1[0:364]);
qnet_interp = np.interp(t/86400, range(364), qnet1[0:364]);
taux_interp = np.interp(t/86400, range(364), taux1[0:364]);
tauy_interp = np.interp(t/86400, range(364), tauy1[0:364]);

dz = .5;
z = np.arange(dz/2,800,dz);

T_shelf =  np.interp(z, -depth_profile, temp_profile);
S_shelf =  np.interp(z, -depth_profile, salt_profile);

t_year = np.arange(0,86400*365,dt)
fwflux_temp = np.interp(t_year/86400, range(364), fwflux1[0:364]);
fwflux_mean = np.mean(fwflux_temp);
fwflux = np.tile(fwflux_temp,int(days/365)+1);
#fwflux = np.concatenate((fwflux_temp,fwflux_temp,fwflux_temp*fwflux_spike,fwflux_temp,fwflux_temp,fwflux_temp,fwflux_temp,fwflux_temp,fwflux_temp),axis = None);


qnet_temp = np.interp(t_year/86400, range(364), qnet1[0:364]);
qnet_mean = np.mean(qnet_temp);
qnet = np.tile(qnet_temp,int(days/365)+1);

taux_temp = np.interp(t_year/86400, range(364), taux1[0:364]);
taux_mean = np.mean(taux_temp);
taux = np.tile(taux_temp,int(days/365)+1);

tauy_temp = np.interp(t_year/86400, range(364), tauy1[0:364]);
tauy_mean = np.mean(tauy_temp);
tauy = np.tile(tauy_temp,int(days/365)+1);

salts = [1];
phi_now = [0.002];

thick_CDW = 800;
thin_CDW = 1598;

h_interface = thick_CDW;
CDWthickness = (1600-h_interface)/2;

T_cdw=1; S_cdw=34.5; T_out=-1; S_out=33.8;

T_start = np.zeros(shape=len(z));
T_start = np.tanh((z-z[h_interface])/100)*(T_cdw-T_out)/2+(T_cdw+T_out)/2;
T_start[h_interface:len(z)]= T_start[h_interface:len(z)]+0.6*((z[h_interface:len(z)]-z[h_interface])/1000);

S_start = np.zeros(shape=len(z));
S_start = np.tanh((z-z[h_interface])/100)*(S_cdw-S_out)/2+(S_cdw+S_out)/2;
S_start[h_interface:len(z)]= S_start[h_interface:len(z)]+0.25*((z[h_interface:len(z)]-z[h_interface])/1000);
now = datetime.now()

days = 0.1*365;
filename1 = '/home/erobo/ThompsonResearch/ModelRuns/' + now.strftime("%Y%m%d_%H%M%S") + '_salt' + str(salts[0]) + '_phi' + str(phi_now[0]) + '_CDWstart' + str(CDWthickness)+'_' + str(days/365) + 'years\r\n';
modelfile=open("ModelsRunning.txt", "a+")
modelfile.write(filename1)


for j in range(len(salts)):
    
    # set up time
    dt = 900; # seconds
    #days = 10*365;
    t = np.arange(0,86400*days,dt)
   
    shape = (int(len(t[0::200])),len(z))
    Ts = np.empty(shape=shape)
    Ss = np.empty(shape=shape)
    #Us = np.empty(shape=shape)
    #Vs = np.empty(shape=shape)
    vel = np.empty(shape=shape)

    shape = (2,len(z))
    Ts_temp = np.empty(shape=shape)
    Ss_temp = np.empty(shape=shape)
    Us_temp = np.empty(shape=shape)
    Vs_temp = np.empty(shape=shape)
    vel_temp = np.empty(shape=shape)

    mld = np.empty(shape=int(len(t[0::200])))
    mld_temp = np.empty(shape=len(z));
    Us_temp[0] = np.zeros(len(z))
    Vs_temp[0] = np.ones(len(z))
    
    #Ts_temp[0] = T_shelf;
    Ts_temp[0] = T_start;

    #Ss_temp[0] = S_shelf;
    Ss_temp[0] = S_start;

    counter = 0;

    start = time.time();
    for i in range(len(t)-1):
            if np.remainder(t[i]/86400,20)==0:
                end = time.time()
                print('Day: ',t[i]/86400,' run time: ',np.around(end-start))

            I = qnet_mean;#qnet[i]; # solar heating
            L = 0; # OLR
            E = 0; # evaporation
            P = fwflux_mean*salts[j];#fwflux[i]; #get_precipitation(t[i]); # precipitation
            tau_x = taux_mean;#taux[i]; # zonal wind stress
            tau_y = tauy_mean;#tauy[i]; # meridional wind stress

            Ts_temp[1],Ss_temp[1],Us_temp[1],Vs_temp[1],mld_temp,vel_temp = PWP.model_timestep(Ts_temp[0],Ss_temp[0],Us_temp[0],Vs_temp[0],z,I,L,E,P,tau_x,tau_y,dt,return_MLD=True,T0=T0,S0=S0,rho0=rho0,alpha=alpha,beta=beta,f=f,advection=True,phi=phi_now,verbose=False,vert_diffusivity=1e-4,return_vel=True)    

            if np.remainder(i,200)==0:
                Ts[counter] = Ts_temp[-1];
                Ss[counter] = Ss_temp[-1];
                #Us[counter] = Us_temp[-1];
                #Vs[counter] = Vs_temp[-1];
                mld[counter] = mld_temp;
                vel[counter] = vel_temp;
                counter += 1;

            Ts_temp[0] = Ts_temp[1];
            Ss_temp[0] = Ss_temp[1];
            Us_temp[0] = Us_temp[1];
            Vs_temp[0] = Vs_temp[1];

    end = time.time()
    print('Final Time: ',np.around(end-start))
    
    #T_start = Ts[-1];
    #S_start = Ss[-1];
    
    filename = '/home/erobo/ThompsonResearch/ModelRuns/' + now.strftime("%Y%m%d_%H%M%S") + '_salt' + str(salts[j]) + '_phi' + str(phi_now[j]) + '_CDWstart' + str(CDWthickness)+'_' + str(days/365) + 'years';
    T_name = filename + '_T.csv';
    S_name = filename + '_S.csv';
    mld_name = filename + '_mld.csv';
    vel_name = filename + '_vel.csv';   
 
    np.savetxt(T_name,Ts[0::2,0::2],fmt='%1.4e',delimiter=',');
    np.savetxt(S_name,Ss[0::2,0::2],fmt='%1.4e',delimiter=',');
    np.savetxt(mld_name,mld[0::2],fmt='%1.4e',delimiter=',');
    np.savetxt(vel_name,vel[0::2,0::2],fmt='%1.4e',delimiter=',');
    print(filename);
