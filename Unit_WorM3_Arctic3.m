clc; clear; 

% Boundary conditions

A = 9.87e11; %m2, Area of the mixed layer in the East Siberian Sea, Jakbosson et al., 2002
ML = 12.8; %m, mixed layer depth of the East Siberian Sea, this study 
V = A*ML; %m3, water volume

% Atmospheric deposition

Dep = 4.1; %g/s, Hg total deposition flux in the Arctic Ocean, Dastoor et al., 2022(SI)
D(1) = 1/3*Dep*10e-6*A/365/86400; % Hg(0) deposition to the East Siberian Sea, Song et al., 2025 2:1 ratio of Hg(II)vsHg(0)
D(2) = 2/3*Dep*10e-6*A/365/86400; % Hg(II) deposition to the East Siberian Sea, Song et al., 2025 2:1 ratio of Hg(II)vsHg(0)


% Meltwater 
ice_b = 3.881; % Sea ice concentration (%) in Sep 1st
ice_a = 3.405; % Sea ice concentration (%) in Sep 15th
ice = (ice_b+ice_a)/200; % Average sea ice concentration
ice_thick = 1.5; % m, Kwok et al., 2018
M_vol = A*(ice_b-ice_a)/100*ice_thick/15/86400; % Meltwater input to the East Siberian Sea, Satellite data
M_THg = 2.6*10e-12*200.59*1000; %g/m3, HgII in sea ice melt, Kim et al., 2023

M_ratio = 0.0142; % Hg0/THg ratio in the Antarctic sea ice, Mastromonaco et al., 2016
M_Hg2 = M_THg*(1-M_ratio);
M_Hg0 = M_THg*M_ratio; 
M(1) = M_vol*M_Hg0;
M(2) = M_vol*M_Hg2;

% River

R_vol_kolyma = 8820; %m3/s, 2024 september, ArcticGRO
R_vol_indigirka = 2652; %m3/s, 1936-1998, september, ArcticcRIMS
R_Hg_kolyma = 23.9*200.59*10e-9; %g/m3, Hg in Kolyma River, Zolko et al., 2020
R_Hg_indigirka = 23.9*200.59*10e-9; %g/m3, 같다고 가정

R_kolyma = R_vol_kolyma*R_Hg_kolyma;
R_indigirka = R_vol_indigirka*R_Hg_indigirka;

R_total = R_kolyma + R_indigirka;
R_ratio = 0.1; %10 - 30% of Hg(0) in THg in river (Leopold et al., 2010)
R(2) = R_total;

% Coastal erosion

Co_ref = 4.36; %ton yr-1, coastal erosion in the East Siberian Sea, Dastoor et al., 2022
Co(2) = Co_ref*10e6/365/86400; %g/s

% Ocean current - Laptev

L_v = 0.0078; %m s-1, current velocity
L_d = 601*1000; %m, latitudinal distance at 140E
L_vol = L_v*ML*L_d; %m3 s-1, volume transport 
L_THg = 0.53*10e-12*200.59*1000; %g/m3, THg, Heimburger et al., 2015
L(2) = L_vol*L_THg;

% Ocean current - CS

CS_v = -0.0041; %m s-1, current velocity
CS_d = 603*1000; %m, latitudinal distance at 140E
CS_vol = CS_v*12.3*CS_d; %m3 s-1, volume transport 
CS_THg = 1.7*10e-12*200.59*1000; %g/m3, THg, Heimburger et al., 2015
CS_ratio = 0.0367; % Hg0/THg ratio
CS_Hg2 = CS_THg*(1-CS_ratio);
CS_Hg0 = CS_THg*CS_ratio; 
CS(1) = CS_vol*CS_Hg0;
CS(2) = CS_vol*CS_Hg2;

% Ocean current - CAO

CAO_v = 0.0114; %m s-1, current velocity
CAO_d = 1205*1000; %m, latitudinal distance at 140E
CAO_t = CAO_v*ML*CAO_d; %m3 s-1, volume transport 

% Air-sea exchange 

Hg_air = 1.05*10^(-9); %g m-3, air Hg(0), this study
Tw = 0.94; % degree C, in situ seawater temperature
H = exp(-2404.3/(Tw+273.15)+6.92); % Dimensionless Henry's law constant
u_10 = 5.01; %m s-1, wind speed at 10 m
p_sw = 1.028-4.97*10^(-5)*Tw-5.575*10^(-6)*Tw^2; %
n_fw = 4.2844*10^(-4)+(0.0157*(Tw+64.993)^(2)-9.13)^(-1); %
n_sw = n_fw*(1.064+6.067*10^-4*Tw-2.753*10^-6*Tw^2); %
pi_w = 2.26;
M_w = 18.01;
V_b = 12.74;
D_sw_Hg = 7.4*10^-10*(273.15+Tw)*sqrt(pi_w*M_w)/n_sw/V_b^0.6;
v_sw = n_sw/p_sw; %
Sc_Hg = v_sw / D_sw_Hg; %Schmidt number of Hg 
k_w = 0.25*u_10^2*(Sc_Hg/660)^(-0.5)*0.01/3600; %m s-1, mass transfer coefficient


%Redox 


PAR = 24.02; %W m-2, intensity
a_443 = 0.1464; %light absorbance at 443 nm
Kd_490 = 0.1616; %diffuse attenuation coefficient at 490 nm
Kd_PAR = 0.00665+0.874*Kd_490-0.00121*Kd_490^(-1); %diffuse attenuation coefficient at PAR

kr_PAR = 0.00133; % h-1 W-1 m2, Ci et al., 2016

kr_PAR_d = kr_PAR*PAR*(1-exp(-Kd_PAR*ML))/Kd_PAR/ML;% depth-averaged kr_PAR
kr_Dark_d = 0.168/24*2^((Tw-25)/10); %depth-averaged kr_dark

RHg_PAR = 0.4852;%Reducible Hg(II) fraction for PAR
k_PAR(2) =kr_PAR_d*RHg_PAR; %Gross Hg(0) production rate (h-1) for PAR


k_Dark(2) = kr_Dark_d;
k_2 =k_PAR(2)+k_Dark(2);
k(2) = k_2/3600;


ko_PAR = 0.0020; % h-1 W-1 m2, Ci et al., 2016
ko_PAR_d = ko_PAR*PAR*(1-exp(-Kd_PAR*ML))/Kd_PAR/ML;% depth-averaged ko_PAR
ko_Dark = 0.081*2^((Tw-25)/10); %depth-averaged ko_dark

k_1 = ko_PAR_d+ko_Dark; %sum
k(1) = k_1/3600;


%Particle - Kim 방법

Kp = 13031668; %Log Kd = 7.12
C_spm = 0.60/10^6;
F_Hg2 = 1-(1/(1+Kp*C_spm));
d_pw = 1500;
d_sw = 1025;
u = 0.00196;
g = 9.8;
r_pw = 5.0*10^(-6);
V_set = 2/9*((d_pw-d_sw)/u)*g*(r_pw)^2 %m s-1, settling velocity


% Water column diffusion - Soerensen 방법 

Dw = 0.8*10^(-4); %eddy diffusion coefficient
Diff_depth_so = 12.8;
Diff_depth_ss = 120;
Diff_d = 20;
Cw(1) = 110*10^(-15)*200.59*1000; % Subsurface Hg(0)
Cw(2) = 1300*10^(-15)*200.59*1000;% Subsurface Hg(II)

% c(1) = aqueous Hg(0); c(2) = aqueous Hg(II)

C0 = [1e-8 1e-7];

options = odeset('RelTol', 1e-13, 'AbsTol', [1e-18 1e-18], 'NonNegative', [1 2]);

tmax =  5*365 * 24 * 60 * 60; % 30 days
tspan = [0 tmax];

[PreIndT, PreIndConc] = ode15s(@Unit_WorM3_ODE_Arctic2,tspan,C0,options,A,V,D,M,R,Co,L,CS,CAO_t,V_set,F_Hg2,k,H,k_w,Hg_air,ice,Cw,Dw,Diff_d);

figure
plot(PreIndT/(5*365*24*60*60), PreIndConc(:,1),'-o')
xlabel('Pre-industrial time (years)')
ylabel('Hg(0) concentration in air (g/m3)')




c_ss = PreIndConc(end,:);

c1_ss = c_ss(1);
c2_ss = c_ss(2);

% ---- Flux 계산 (g s^-1) ----

Diff_Hg0      = Dw*(Cw(1) - c1_ss)/Diff_d*A;
Diff_HgII      = Dw*(Cw(2) - c2_ss)/Diff_d*A;
AS = (1-ice)*k_w*(c1_ss-Hg_air/H)*A;
CAO_Hg0       = CAO_t * c1_ss;
CAO_HgII       = CAO_t * c2_ss;
S0_HgIISink    = V_set*A*F_Hg2 * c2_ss;
S0_Hg0toHgII   = k(1)*V * c1_ss;
S0_HgIItoHg0   = k(2)*V * c2_ss;
D_HgII = D(2)*(1-ice);




