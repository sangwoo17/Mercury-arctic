clc; clear; 

% Boundary conditions

A = 6.71154e11; %m2, Area of the mixed layer in the CAO, rstudio
ML = 14.7725*1.001; %m, mixed layer depth of the Makarov Basin
V = A*ML; %m3, water volume
ice = 1; % sea ice concentration                                           %% 변경 

% Atmospheric deposition

Dep = 4.1; %ug m-2 y-1, Fall, Hg total deposition flux in the Arctic Ocean, Dastoor et al., 2022(SI)
D(2) = (1-ice)*Dep*1e-6*A/365/86400; %g s-1, Hg(II) deposition

% Meltwater 
ice_b = 0.4348; % Sea ice concentration (%) in Sep 1st
ice_a = 0.6489; % Sea ice concentration (%) in Sep 30th
ice_thick = 1.5; % m, Kwok et al., 2018
M_vol = A*(ice_b-ice_a)*ice_thick/30/86400; % Meltwater input to the East Siberian Sea, Satellite data

% Ocean current - Makarov Basin 

MB_v = 0.0299; %m s-1, current velocity
MB_d = 222.639*1e3; %m, latitudinal distance at 140E
MB_vol = MB_v*ML*MB_d; %m3 s-1, volume transport 
MB_THg = 1.6*1e-12*200.59*1e3; %g/m3, THg in the Makarov Basin, Eom et al., 2025
MB_ratio = 0.07535; % Hg0/THg ratio from the ESS and CS in this study and Kim et al., 2020
MB_Hg2 = MB_THg*(1-MB_ratio);
MB_Hg0 = MB_THg*MB_ratio; 
MB(1) = MB_vol*MB_Hg0;
MB(2) = MB_vol*MB_Hg2;

% Ocean current - CS

CS_v = 0.0273; %m s-1, current velocity
CS_d = 1104*1000; %m, latitudinal distance at 140E
CS_vol = CS_v*12.3*CS_d; %m3 s-1, volume transport 
CS_THg = 1.7*1e-12*200.59*1e3; %g/m3, THg, Kim et al., 2020
CS_ratio = 0.02647; % Hg0/THg ratio
CS_Hg2 = CS_THg*(1-CS_ratio);
CS_Hg0 = CS_THg*CS_ratio; 
CS(1) = CS_vol*CS_Hg0;
CS(2) = CS_vol*CS_Hg2;

% Ocean current - ESS

ESS_v = 0.0111; %m s-1, current velocity
ESS_d = 1256.5*1000; %m, latitudinal distance at 140E
ESS_vol = ESS_v*12.8*ESS_d; %m3 s-1, volume transport 
ESS_THg = 1.0*1e-12*200.59*1e3; %g/m3, THg, Kim et al., 2020
ESS_ratio = 0.04834; % Hg0/THg ratio
ESS_Hg2 = ESS_THg*(1-ESS_ratio);
ESS_Hg0 = ESS_THg*ESS_ratio; 
ESS(1) = ESS_vol*ESS_Hg0;
ESS(2) = ESS_vol*ESS_Hg2;

% Ocean current - CAO - North

CAO_nv = 0.0174; %m s-1, current velocity
CAO_nd = 1051.26*1000; %m, latitudinal distance at 140E
CAO_nt = CAO_nv*ML*CAO_nd; %m3 s-1, volume transport 

% Ocean current - CAO - East

CAO_ev = 0.0048*1.001; %m s-1, current velocity
CAO_ed = 667.9169*1000; %m, latitudinal distance at 140E
CAO_et = CAO_ev*ML*CAO_ed; %m3 s-1, volume transport 

% Air-sea exchange 

Hg_air = 1.11*1e-9; %g m-3, air Hg(0), this study                                     %% 변경
Tw = -0.35; % degree C, in situ seawater temperature                                  %% 변경
H = exp(-2404.3/(Tw+273.15)+6.92); % Dimensionless Henry's law constant
u_10 = 7.39; %m s-1, wind speed at 10 m                                               %% 변경
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

%Particle - Kim 방법

Kp = 10^(6.83); %
C_spm = 0.02/10^6;
F_Hg2 = 1-(1/(1+Kp*C_spm));
d_pw = 1500;
g = 9.8;
r_pw = 5*10^(-6)*1.001;
V_set = 2/9*((d_pw-p_sw*1000)/(n_sw*0.1))*g*(r_pw)^2; %m s-1, settling velocity

%Redox 

UV_B = 0.02*(1-ice); %W m-2, intensity
UV_A = 4.2*(1-ice) ; %W m-2, intensity 
PAR = 24.02*(1-ice); %W m-2, intensity
a_443 = 0.1464; %light absorbance at 443 nm
Kd_490 = 0.1616; %diffuse attenuation coefficient at 490 nm
Kd_380 = 3.313*(a_443)^2+0.887*a_443+0.03; %diffuse attenuation coefficient at 380 nm
Kd_340 = -1.965*(Kd_380)^2+1.979*Kd_380-0.009; %diffuse attenuation coefficient at 340 nm
Kd_305 = -7.07*(Kd_380)^2+4.424*Kd_380-0.008; %diffuse attenuation coefficient at 305 nm
Kd_PAR = 0.00665+0.874*Kd_490-0.00121*Kd_490^(-1); %diffuse attenuation coefficient at PAR

UV_B_exp = 1.4; % W m-2, UV-B intensity in this study
UV_A_exp = 6.2; % W m-2, UV-A intensity in this study
kr_UVA = 0.0806/UV_A_exp; %h-1 W-1 m2, light-normalized kr_UVA              %% 변경
kr_UVB = 0.208/UV_B_exp; %h-1 W-1 m2, light-normalized kr_UVB              %% 변경 

kr_UVB_d = kr_UVB*UV_B*(1-exp(-Kd_305*ML))/Kd_305/ML; % depth-averaged kr_UVB
kr_UVA_d = kr_UVA*UV_A*(1-exp(-Kd_340*ML))/Kd_340/ML; % depth-averaged kr_UVA
kr_Dark_d = 0.168/24*2^((Tw-25)/10); %depth-averaged kr_dark
RHg_UVA = 0.525*1/(1+Kp*C_spm); %Reducible Hg(II) fraction for UVA;                       %% 변경
RHg_UVB = 0.8842*1/(1+Kp*C_spm); %Reducible Hg(II) fraction for UVB;                       %% 변경

k_UVB(2) = kr_UVB_d*RHg_UVB; %Gross Hg(0) production rate (h-1) for UVB
k_UVA(2) = kr_UVA_d*RHg_UVA; %Gross Hg(0) production rate (h-1) for UVA

k_Dark(2) = kr_Dark_d;
k_2 = k_UVB(2)+k_UVA(2)+k_Dark(2);
k(2) = k_2/3600;

UV_B_exp2 = 1.35;
UV_A_exp2 = 5.33;
ko_UVA = 0.3852/UV_A_exp2; %h-1 W-1 m2, light-normalized kr_UVA              %% 변경
ko_UVB = 1.909/UV_B_exp2; %h-1 W-1 m2, light-normalized kr_UVB               %% 변경
ko_UVB_d = ko_UVB*UV_B*(1-exp(-Kd_305*ML))/Kd_305/ML; % depth-averaged kr_UVB
ko_UVA_d = ko_UVA*UV_A*(1-exp(-Kd_340*ML))/Kd_340/ML; % depth-averaged kr_UVA
k_UVB(1) = ko_UVB_d; %Gross Hg(0) production rate (h-1) for UVB
k_UVA(1) = ko_UVA_d; %Gross Hg(0) production rate (h-1) for UVA
ko_Dark = 0.068*2^((Tw-25)/10); %depth-averaged ko_dark

k_1 = k_UVB(1)+k_UVA(1)+ko_Dark; %sum
k(1) = k_1/3600;


% Water column diffusion - Soerensen 방법 

Dw = 0.8*10^(-4); %eddy diffusion coefficient
Diff_d = 100/2;
Cw(1) = 252.6*10^(-15)*200.59*1000; % Subsurface Hg(0)
Cw(2) = 996.7*10^(-15)*200.59*1000;% Subsurface Hg(II)

% c(1) = aqueous Hg(0); c(2) = aqueous Hg(II)

C0 = [1e-8 1e-7];

options = odeset('RelTol', 1e-13, 'AbsTol', [1e-18 1e-18], 'NonNegative', [1 2]);

tmax =  2*365 * 24 * 60 * 60; % 30 days
tspan = [0 tmax];

[PreIndT, PreIndConc] = ode15s(@Unit_WorM3_ODE_CAO,tspan,C0,options,A,V,D,M_vol,MB,CS,ESS,CAO_nt,CAO_et,V_set,F_Hg2,k,H,k_w,Hg_air,ice,Cw,Dw,Diff_d);

%figure
%plot(PreIndT/(365*24*60*60), PreIndConc(:,1),'-o')
%xlabel('Pre-industrial time (years)')
%ylabel('Hg(0) concentration in air (g/m3)')

c_ss = PreIndConc(end,:);

c1_ss = c_ss(1);
c2_ss = c_ss(2);

% ---- Flux 계산 (g s^-1) ----

a_Diff_Hg0      = Dw*(Cw(1) - c1_ss)/Diff_d*A;
a_Diff_HgII      = Dw*(Cw(2) - c2_ss)/Diff_d*A;
a_AS = (1-ice)*k_w*(c1_ss-Hg_air/H)*A;
a_CAOn_Hg0       = CAO_nt * c1_ss;
a_CAOn_HgII       = CAO_nt * c2_ss;
a_CAOe_Hg0       = CAO_et * c1_ss;
a_CAOe_HgII       = CAO_et * c2_ss;
a_S0_HgIISink    = V_set*A*F_Hg2 * c2_ss;
a_S0_Hg0toHgII   = k(1)*V * c1_ss;
a_S0_HgIItoHg0   = k(2)*V * c2_ss;
a_D_HgII = D(2)*(1-ice);
a_M_HgII = M_vol*c2_ss;
a_M_Hg0 = M_vol*c1_ss;
a_MB(1) = MB(1);
a_MB(2) = MB(2);
a_CS(1) = CS(1);
a_CS(2) = CS(2);
a_ESS(1) = ESS(1);
a_ESS(2) = ESS(2);

FluxName = {
    'Diff_Hg0'
    'Diff_HgII'
    'AirSea_Hg0'
    'CAOn_Hg0'
    'CAOn_HgII'
    'CAOe_Hg0'
    'CAOe_HgII'
    'Settling_HgII'
    'Hg0_to_HgII'
    'HgII_to_Hg0'
    'Deposition_HgII'
    'Melt_HgII'
    'Melt_Hg0'
    'MB_Hg0'
    'MB_HgII'
    'CS_Hg0'
    'CS_HgII'
    'ESS_Hg0'
    'ESS_HgII'
    'c(1)'
    'c(2)'
};

FluxValue = [
    a_Diff_Hg0
    a_Diff_HgII
    a_AS
    a_CAOn_Hg0
    a_CAOn_HgII
    a_CAOe_Hg0
    a_CAOe_HgII
    a_S0_HgIISink
    a_S0_Hg0toHgII
    a_S0_HgIItoHg0
    a_D_HgII
    a_M_HgII
    a_M_Hg0
    a_MB(1)
    a_MB(2)
    a_CS(1)
    a_CS(2)
    a_ESS(1)
    a_ESS(2)
    c1_ss
    c2_ss
];

T = table(FluxName, FluxValue);


