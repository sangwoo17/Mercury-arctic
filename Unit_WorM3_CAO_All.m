clc; clear; 

sp_values = 27.4:0.1:29.1;

ResultTable = table();

for s_p = sp_values

% Boundary conditions

s_bg = 34.8;
% s_p = 31;
Tw = 0.1261*(s_p)^2-7.7555*(s_p)+117.57; % degree C, in situ seawater temperature, ITS-90 기준                      %% 변경
t_68 = (Tw+273.15-0.00025*273.15)/(1-0.00025)-273.15; % degree C, in situ seawater temperature, IPTS-68 기준
p_sw_A = 0.824493 - 4.0899*10^(-3) * t_68 + 7.6438 * 10^(-5) * (t_68)^2 -8.2467 * 10^(-7) *(t_68)^3 + 5.3875*10^(-9) * (t_68)^4;
p_sw_B = -5.72466*10^(-3)+1.0227*10^(-4)*t_68-1.6546*10^(-6)*(t_68)^2;
p_sw_C = 4.8314*10^(-4);
p_w = 999.842594 + 6.793952*10^(-2)*t_68 - 9.09529*10^(-3)*(t_68)^2 + 1.001685*10^(-4)*(t_68)^3 - 1.120083*10^(-6)*(t_68)^4 + 6.536336*10^(-9)*(t_68)^5;
p_sw = (p_w +  p_sw_A * s_p + p_sw_B * s_p^(3/2) + p_sw_C * s_p)/1000;
p_uhw = 1024.8; % kg m-3, density of UHW in the Makarov Basin 
MLD = 13 * (p_uhw - p_sw*1000)^(-0.45); %m, mixed layer depth of the Makarov Basin
f_ow = s_p / s_bg;
f_meteo = 1-f_ow;

A = 6.891e11; %m2, Area of the mixed layer in the CAO, rstudio
V = A*MLD; %m3, water volume
ice = 0; % sea ice concentration                                           %% 변경 0
ABL = 236; % Atmospheric boundary layer, Peng et al., 2023, Near-neutral ABL (NBL)
V_A = A*ABL; % Atmospheric volume
air_ratio = 0.00733; % ratio of HgII to Hg0, 0.00733

% Atmospheric transport - East Siberian Sea
L_ess_wind = 3.1475; % m s-1, wind speed 
L_ess_dist = 1038.32*1000; % m, longitudinal distance
L_ess_v = L_ess_wind * L_ess_dist * ABL; % m3 s-1, volume transport from the ESS
air_Hg_ess = 1.05; % TGM in the ESS
L_ess(1) = L_ess_v*air_Hg_ess*10^-9; % g s-1, Volume transport of Hg0

% Atmospheric transport - Chukchi Sea
L_cs_wind = -2.5303; % m s-1, wind speed
L_cs_dist = 610.01 * 1000; %m, longitudinal distance
L_cs_v = L_cs_wind * L_cs_dist * ABL; % m3 s-1, volume transport from the CS

% Atmospheric transport - West 
L_w_wind = -4.5416; % m s-1, wind speed
L_w_dist = 500 * 1000; %m, longitudinal distance
L_w_v =  L_w_wind * L_w_dist * ABL; % m3 s-1, volume transport from the west

% Atmospheric transport - North 

L_n_wind = -2.4825; % m s-1, wind speed
L_n_dist = 1105 * 1000; %m, longitudinal distance
L_n_v =  L_n_wind * L_n_dist * ABL; % m3 s-1, volume transport from the west
air_Hg_n = 1.176829; % TGM in the northern CAO
L_n(1) = L_n_v*air_Hg_n*10^-9; % g s-1, Volume transport of Hg0

% Atmospheric transport - East
L_e_wind = -2.5448; % m s-1, wind speed
L_e_dist = 500 * 1000; %m, longitudinal distance
L_e_v =  L_e_wind * L_e_dist * ABL; % m3 s-1, volume transport from the west
air_Hg_e = 1.144374; % TGM in the CB 
L_e(1) = L_e_v*air_Hg_e*10^-9; % g s-1, Volume transport of Hg0

% Atmospheric transport - vertical
L_v_v = L_ess_v + L_cs_v + L_w_v - L_e_v - L_n_v;

% Dry deposition
v_dep = 1.5; %cm s-1, Orbist et al., 2017
V_dep = v_dep*0.01*A; % m3 s-1, volume flux 

% Precipitation
precip = 0.001061/24/3600*A; % m3 s-1, volume transport of precipitation 
precip_thg = 3.1*10^-12*200.59*1000; %g m-3, rainwater THg
P = (1-ice)*precip * precip_thg; % g s-1, wet deposition

% Meltwater 
ice_b = 176.07*10^9; % Sea ice volume (m3) in Sep 1st
ice_a = 278.32*10^9; % Sea ice volume (m3) in Sep 30th
M_vol = (ice_b-ice_a)/30/86400; % Meltwater input to the East Siberian Sea, Satellite data

% Ocean current - Makarov Basin 
MB_v = 0.022997; %m s-1, current velocity
MB_d = 500377.2; %m, latitudinal distance at 140E
MB_vol = MB_v*MLD*MB_d; %m3 s-1, volume transport 
MB_THg = 1.6*1e-12*200.59*1e3; %g/m3, THg in the Makarov Basin, Eom et al., 2025
MB_ratio = 0.07535; % Hg0/THg ratio from the ESS and CS in this study and Kim et al., 2020
MB_Hg2 = MB_THg*(1-MB_ratio);
MB_Hg0 = MB_THg*MB_ratio; 
MB(1) = MB_vol*MB_Hg0*(1-f_meteo);
MB(2) = MB_vol*MB_Hg2*(1-f_meteo);

% Ocean current - CS
CS_v = 0.035393; %m s-1, current velocity
CS_d = 610011.6; %m, latitudinal distance at 140E
CS_vol = CS_v*MLD*CS_d; %m3 s-1, volume transport 
CS_THg = 1.7*1e-12*200.59*1e3; %g/m3, THg, Kim et al., 2020
CS_ratio = 0.02647; % Hg0/THg ratio
CS_Hg2 = CS_THg*(1-CS_ratio);
CS_Hg0 = CS_THg*CS_ratio; 
CS(1) = CS_vol*CS_Hg0*(1-f_meteo);
CS(2) = CS_vol*CS_Hg2*(1-f_meteo);

% Ocean current - ESS
ESS_v = 0.011518; %m s-1, current velocity
ESS_d = 1038317.6; %m, latitudinal distance at 140E
ESS_vol = ESS_v*MLD*ESS_d; %m3 s-1, volume transport 
ESS_THg = 1.0*1e-12*200.59*1e3; %g/m3, THg, Kim et al., 2020
ESS_ratio = 0.04834; % Hg0/THg ratio
ESS_Hg2 = ESS_THg*(1-ESS_ratio);
ESS_Hg0 = ESS_THg*ESS_ratio; 
ESS(1) = ESS_vol*ESS_Hg0*(1-f_meteo);
ESS(2) = ESS_vol*ESS_Hg2*(1-f_meteo);

% Ocean current - CAO - North
CAO_nv = 0.016198; %m s-1, current velocity
CAO_nd = 1104564.6; %m, latitudinal distance at 140E
CAO_nt = CAO_nv*MLD*CAO_nd; %m3 s-1, volume transport 

% Ocean current - CAO - East
CAO_ev = 0.029777; %m s-1, current velocity
CAO_ed = 500377.2; %m, latitudinal distance at 140E
CAO_et = CAO_ev*MLD*CAO_ed; %m3 s-1, volume transport 

% Freshwater
fr_vol = (ESS_vol+CS_vol+MB_vol)*f_meteo;
fr_THg = 5.46*1e-12*200.59*1e3; %g/m3, THg in the Makarov Basin, Eom et al., 2025       5.46
fr_ratio = 0.04834; % median of the CS, MB, and ESS;
fr_Hg2 = fr_THg*(1-fr_ratio);
fr_Hg0 = fr_THg*fr_ratio;
fr(1) = fr_vol*fr_Hg0;
fr(2) = fr_vol*fr_Hg2;

% Air-sea exchange                       
H = exp(-2404.3/(Tw+273.15)+6.92); % Dimensionless Henry's law constant
u_10 = 7.39; %m s-1, wind speed at 10 m                                               %% 변경
n_fw = 4.2844*10^(-5)+(0.157*(Tw+64.993)^(2)-91.296)^(-1); % kg m-1 s-1, freshwater dynamic viscosity 
n_sw_A = 1.541 + 1.998*10^(-2)*Tw-9.52*10^(-5)*Tw^2;
n_sw_B = 7.974 - 7.561*10^(-2)*Tw+4.724*10^(-4)*Tw^2;
n_sw = 10*n_fw*(1+n_sw_A*(s_p)/1000+n_sw_B*(s_p/1000)^2);% g cm s-1, seawater dynamic viscosity 
pi_w = 2.26;
M_w = 18.01;
V_b = 12.74;
D_sw_Hg = 7.4*10^-10*(273.15+Tw)*sqrt(pi_w*M_w)/n_sw/V_b^0.6;
v_sw = n_sw/p_sw; %
Sc_Hg = v_sw / D_sw_Hg; %Schmidt number of Hg 
k_w = 0.25*u_10^2*(Sc_Hg/660)^(-0.5)*0.01/3600; %m s-1, mass transfer coefficient

%Particle - Kim 방법
Kp = 10^(6.5); % Cui et al., 2021
C_spm = 21.4*10^-9; % Phoebe J. Lam, Size-fractionated particle composition, Arctic Ocean
F_Hg2 = 1-(1/(1+Kp*C_spm));
V_set = 24.9*10^(-3)/24/3600*A; % g s-1, particle settling velocity, fall averaged without outliers

%Redox 

UV_B_in = 0.0061;
UV_A_in = 1.02;
PAR_in = 16.2;
I_trans = 0.028; % Median light transmittance, Katlein et al., 2019

UV_B = UV_B_in*(1-ice) + UV_B_in*ice*I_trans; %W m-2, intensity
UV_A = UV_A_in*(1-ice) + UV_A_in*ice*I_trans; %W m-2, intensity 
PAR = PAR_in*(1-ice) + PAR_in*ice*I_trans; ; %W m-2, intensity
a_443 = 0.06; %light absorbance at 443 nm
Kd_490 = 0.08; %diffuse attenuation coefficient at 490 nm
Kd_380 = 3.313*(a_443)^2+0.887*a_443+0.03; %diffuse attenuation coefficient at 380 nm
Kd_340 = -1.965*(Kd_380)^2+1.979*Kd_380-0.009; %diffuse attenuation coefficient at 340 nm
Kd_305 = -7.07*(Kd_380)^2+4.424*Kd_380-0.008; %diffuse attenuation coefficient at 305 nm
Kd_PAR = 0.00665+0.874*Kd_490-0.00121*Kd_490^(-1); %diffuse attenuation coefficient at PAR

%RHg_UVA = 0.522*1/(1+Kp*C_spm); %Reducible Hg(II) fraction for UVA;
RHg_UVA = (-0.0196*(s_p)^2+1.1782*(s_p)-17.153)*1/(1+Kp*C_spm); %Reducible Hg(II) fraction for UVA;                      %% 변경
RHg_UVB = 0.854*1/(1+Kp*C_spm); %Reducible Hg(II) fraction for UVB;                       %% 변경
RHg_PAR = 0.4852*1/(1+Kp*C_spm); %Reducible Hg(II) fraction for PAR;   


UV_B_exp = 1.4; % W m-2, UV-B intensity in this study
UV_A_exp = 6.2; % W m-2, UV-A intensity in this study
%kr_UVA = 0.081/UV_A_exp; %h-1 W-1 m2, light-normalized kr_UVA              %% 변경
kr_UVA = (0.0096*(s_p)^2 - 0.5586*(s_p)+8.2102)/UV_A_exp;
kr_UVB = 0.199/UV_B_exp;
%kr_UVB = (0.8495*(RHg_UVA)^2-0.4843*(RHg_UVA)+0.2377)/UV_B_exp; %h-1 W-1 m2, light-normalized kr_UVB              %% 변경 
kr_PAR = 0.000934;
%kr_PAR = 0.00133; %h-1 W-1 m2, light-normalized kr_PAR              %% 변경 

kr_UVB_d = kr_UVB*UV_B*(1-exp(-Kd_305*MLD))/Kd_305/MLD; % depth-averaged kr_UVB
kr_UVA_d = kr_UVA*UV_A*(1-exp(-Kd_340*MLD))/Kd_340/MLD; % depth-averaged kr_UVA
kr_PAR_d = kr_PAR*PAR*(1-exp(-Kd_PAR*MLD))/Kd_PAR/MLD; % depth-averaged kr_UVA
kr_Dark_d = 0.168/24*(2^((Tw-25)/10)); %depth-averaged kr_dark

k_UVB(2) = kr_UVB_d*RHg_UVB; %Gross Hg(0) production rate (h-1) for UVB
k_UVA(2) = kr_UVA_d*RHg_UVA; %Gross Hg(0) production rate (h-1) for UVA
k_PAR(2) = kr_PAR_d*RHg_PAR;

k_Dark(2) = kr_Dark_d;
k_2 = k_UVB(2)+k_UVA(2)+k_Dark(2);%+k_PAR(2);

k(2) = k_2/3600;

UV_B_exp2 = 1.35;
UV_A_exp2 = 5.33;
%ko_UVA = (0.0022*Tw^2-0.0321*Tw+0.3455)/UV_A_exp2; %h-1 W-1 m2, light-normalized kr_UVA              %% 변경
ko_UVA = 0.38/UV_A_exp2;
ko_UVB = 1.9/UV_B_exp2; %h-1 W-1 m2, light-normalized kr_UVB               %% 변경
ko_PAR = 0.002; %h-1 W-1 m2, light-normalized ko_PAR               %% 변경

ko_UVB_d = ko_UVB*UV_B*(1-exp(-Kd_305*MLD))/Kd_305/MLD; % depth-averaged kr_UVB
ko_UVA_d = ko_UVA*UV_A*(1-exp(-Kd_340*MLD))/Kd_340/MLD; % depth-averaged kr_UVA
ko_PAR_d = ko_PAR*PAR*(1-exp(-Kd_PAR*MLD))/Kd_PAR/MLD; % depth-averaged kr_UVA
k_UVB(1) = ko_UVB_d; %Gross Hg(0) production rate (h-1) for UVB
k_UVA(1) = ko_UVA_d; %Gross Hg(0) production rate (h-1) for UVA
k_PAR(1) = ko_PAR_d; %Gross Hg(0) production rate (h-1) for PAR
ko_Dark = 0.068*(2^((Tw-25)/10)); %depth-averaged ko_dark

k_1 = k_UVB(1)+k_UVA(1)+ko_Dark; %sum
k(1) = k_1/3600;

% Water column diffusion - Soerensen 방법 
Diff_d = 35/2;
g = 9.8;
r = 0.2; % empirical coefficient, Randelhoff et al., 
e_diss = 1 * 10^(-8); % W kg-1, dissipation rate of turbulent kinetic energy, SML
N_bv =  g / (p_sw*1000) * (p_uhw - p_sw*1000) / Diff_d; % s-2, Brunt-Vaisala frequency
Dw = r * e_diss * (N_bv)^(-1);  %eddy diffusion coefficient

Cw(1) = 298.462 *10^(-15)*200.59*1000; % Subsurface Hg(0)
Cw(2) = 726.154*10^(-15)*200.59*1000;% Subsurface Hg(II)

% c(1) = aqueous Hg(0); c(2) = aqueous Hg(II) c(3) = atmospheric Hg(0); c(4) = atmospheric Hg(II)
C0 = [1e-8 1e-7 1e-9];

options = odeset('RelTol', 1e-13, 'AbsTol', [1e-18 1e-18 1e-18], 'NonNegative', [1 2 3]);

tmax =  365 * 24 * 60 * 60; % 30 days
tspan = [0 tmax];

[PreIndT, PreIndConc] = ode15s(@Unit_WorM3_ODE_CAO_All,tspan,C0,options,A,V,M_vol,MB,CS,ESS,CAO_nt,CAO_et,V_set,F_Hg2,Kp,k,H,k_w,ice,Cw,Dw,Diff_d,V_A,L_ess,L_cs_v,L_w_v,L_n,L_e,L_v_v,V_dep,air_ratio,P,fr);

%figure
%plot(PreIndT/(365*24*60*60), PreIndConc(:,1),'-o')
%xlabel('Time (years)')
%ylabel('Hg(0) in water (g/m3)')

%figure
%plot(PreIndT/(365*24*60*60), PreIndConc(:,2),'-o')
%xlabel('Time (years)')
%ylabel('Hg(II) in water (g/m3)')

%figure
%plot(PreIndT/(365*24*60*60), PreIndConc(:,3),'-o')
%xlabel('Time (years)')
%ylabel('Hg(0) in air (g/m3)')

c_ss = PreIndConc(end,:);

c1_ss = c_ss(1);
c2_ss = c_ss(2);
c3_ss = c_ss(3);

% ---- Flux 계산 (g s^-1) ----

a_MLD = MLD;
a_Dw = Dw;
a_Diff_Hg0      = Dw*(Cw(1) - c1_ss)/Diff_d*A;
a_Diff_HgII      = Dw*(Cw(2) - c2_ss)/Diff_d*A;
a_AS               = (1-ice)*k_w*(c1_ss-c3_ss/H)*A;
a_Dep               = (1-ice)* V_dep * c3_ss*air_ratio;
a_CAOn(1)       = CAO_nt * c1_ss;
a_CAOn(2)       = CAO_nt * c2_ss;
a_CAOe(1)       = CAO_et * c1_ss;
a_CAOe(2)       = CAO_et * c2_ss;
a_S0_HgIISink    = V_set*Kp*10^(-6)* c2_ss;
a_S0_Hg0toHgII   = k(1)*V * c1_ss;
a_S0_HgIItoHg0   = k(2)*V * c2_ss;
a_M(2) = M_vol*c2_ss;
a_M(1) = M_vol*c1_ss;
a_MB(1) = MB(1);
a_MB(2) = MB(2);
a_CS(1) = CS(1);
a_CS(2) = CS(2);
a_ESS(1) = ESS(1);
a_ESS(2) = ESS(2);
%a_L_ess(1) = L_ess(1);
%a_L_n(1) = L_n(1);
%a_L_e(1) = L_e(1);
%a_L_w_v(1) = L_w_v*c3_ss;
%a_L_cs_v(1) = L_cs_v*c3_ss;
%a_L_v(1) = -L_v_v * c3_ss;
%a_L(1) = (L_ess(1) - L_n(1) - L_e(1) + L_w_v*c3_ss + L_cs_v*c3_ss -L_v_v*c3_ss);
a_P = P;
a_fr(1) = fr_vol*fr_Hg0;
a_fr(2) = fr_vol*fr_Hg2;


NewRow = table( ...
    s_p, c1_ss, c2_ss, c3_ss, MLD, Dw, ...
    a_Diff_Hg0, a_Diff_HgII, a_AS, a_Dep, ...
    a_CAOn(1), a_CAOn(2), a_CAOe(1), a_CAOe(2), ...
    a_S0_HgIISink, a_S0_Hg0toHgII, a_S0_HgIItoHg0, ...
    a_M(2), a_M(1), a_MB(1), a_MB(2), a_CS(1), a_CS(2), a_ESS(1), a_ESS(2), a_P, a_fr(1), a_fr(2), ...
    'VariableNames', { ...
    'Salinity', ...
    'Hg0_water', ...
    'HgII_water', ...
    'Hg0_air', ...
    'MixedLayerDepth', ...
    'EddyDiffusivity', ...
    'Diffusion_Hg0', ...
    'Diffusion_HgII', ...
    'AirSeaFlux', ...
    'DryDeposition', ...
    'CAOn_Hg0', ...
    'CAOn_HgII', ...
    'CAOe_Hg0', ...
    'CAOe_HgII', ...
    'P_sink', ...
    'Ox', ...
    'Red', ...
    'M_HgII', ...
    'M_Hg0', ...
    'MB_Hg0', ...
    'MB_HgII', ...
    'CS_Hg0', ...
    'CS_HgII', ...
    'ESS_Hg0', ...
    'ESS_HgII', ...
    'Precip', ...
    'Fr_Hg0', ...
    'Fr_HgII' ...
    } ...
);

ResultTable = [ResultTable; NewRow];

end


writetable(ResultTable,'outputs/Hg_budget_summary.csv')