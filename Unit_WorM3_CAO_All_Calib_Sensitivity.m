clc; clear;

data = readtable('inputs/DGM_Eval_CAO.csv');

eps_rel = 0.001;   % 0.1%

SensTable = table();

for i = 1:height(data)

    s_p_raw = data.Salinity(i);

    if s_p_raw <= 27.5
        continue
    end

    % baseline inputs
    s_p  = s_p_raw;
    Tw   = data.Temperature(i);   % deg C
    ice  = data.ICE(i);           % 0-1
    u_10 = data.WS(i);
    lat  = data.Lat(i);

    GEM_obs = data.GEM(i);
    DGM_obs = data.DGM(i);

    % -------------------
    % baseline run
    % -------------------
    out0 = run_CAO_case(s_p, Tw, ice, u_10, lat, GEM_obs, DGM_obs);

    Y0_DGM  = out0.Hg0_water;
    Y0_AS   = out0.AS;
    Y0_Net  = out0.Net;

    % =========================
    % Temperature perturbation
    % % on Kelvin, not Celsius
    % =========================
    Tk0   = Tw + 273.15;
    Tw_p  = Tk0*(1+eps_rel) - 273.15;
    Tw_m  = Tk0*(1-eps_rel) - 273.15;

    outTp = run_CAO_case(s_p, Tw_p, ice, u_10, lat, GEM_obs, DGM_obs);
    outTm = run_CAO_case(s_p, Tw_m, ice, u_10, lat, GEM_obs, DGM_obs);

    % =========================
    % Salinity perturbation
    % =========================
    s_p_p = s_p*(1+eps_rel);
    s_p_m = s_p*(1-eps_rel);

    if s_p_m <= 26.46
        continue
    end

    outSp = run_CAO_case(s_p_p, Tw, ice, u_10, lat, GEM_obs, DGM_obs);
    outSm = run_CAO_case(s_p_m, Tw, ice, u_10, lat, GEM_obs, DGM_obs);

    % =========================
    % Sea ice perturbation
    % =========================
    ice_p = min(1, ice*(1+eps_rel));
    ice_m = max(0, ice*(1-eps_rel));

    outIp = run_CAO_case(s_p, Tw, ice_p, u_10, lat, GEM_obs, DGM_obs);
    outIm = run_CAO_case(s_p, Tw, ice_m, u_10, lat, GEM_obs, DGM_obs);

    % =========================
    % Wind perturbation
    % =========================
    u_p = u_10*(1+eps_rel);
    u_m = u_10*(1-eps_rel);

    outWp = run_CAO_case(s_p, Tw, ice, u_p, lat, GEM_obs, DGM_obs);
    outWm = run_CAO_case(s_p, Tw, ice, u_m, lat, GEM_obs, DGM_obs);

    % =========================
    % Sensitivity coefficients
    % =========================
    S_T_DGM = ((outTp.Hg0_water - outTm.Hg0_water) / Y0_DGM) / (2*eps_rel);
    S_S_DGM = ((outSp.Hg0_water - outSm.Hg0_water) / Y0_DGM) / (2*eps_rel);
    S_I_DGM = ((outIp.Hg0_water - outIm.Hg0_water) / Y0_DGM) / (2*eps_rel);
    S_W_DGM = ((outWp.Hg0_water - outWm.Hg0_water) / Y0_DGM) / (2*eps_rel);

    S_T_AS  = ((outTp.AS - outTm.AS) / Y0_AS) / (2*eps_rel);
    S_S_AS  = ((outSp.AS - outSm.AS) / Y0_AS) / (2*eps_rel);
    S_I_AS  = ((outIp.AS - outIm.AS) / Y0_AS) / (2*eps_rel);
    S_W_AS  = ((outWp.AS - outWm.AS) / Y0_AS) / (2*eps_rel);

    S_T_Net = ((outTp.Net - outTm.Net) / Y0_Net) / (2*eps_rel);
    S_S_Net = ((outSp.Net - outSm.Net) / Y0_Net) / (2*eps_rel);
    S_I_Net = ((outIp.Net - outIm.Net) / Y0_Net) / (2*eps_rel);
    S_W_Net = ((outWp.Net - outWm.Net) / Y0_Net) / (2*eps_rel);

    NewRow = table(i, s_p, Tw, ice, u_10, ...
        Y0_DGM, Y0_AS, Y0_Net, ...
        S_T_DGM, S_S_DGM, S_I_DGM, S_W_DGM, ...
        S_T_AS,  S_S_AS,  S_I_AS,  S_W_AS, ...
        S_T_Net, S_S_Net, S_I_Net, S_W_Net, ...
        'VariableNames', { ...
        'RowID','Salinity','Temperature','SeaIce','WindSpeed', ...
        'DGM_base','AS_base','Net_base', ...
        'Sens_T_DGM','Sens_S_DGM','Sens_Ice_DGM','Sens_W_DGM', ...
        'Sens_T_AS','Sens_S_AS','Sens_Ice_AS','Sens_W_AS', ...
        'Sens_T_Net','Sens_S_Net','Sens_Ice_Net','Sens_W_Net'});

    SensTable = [SensTable; NewRow];
end

writetable(SensTable, 'outputs/Hg_sensitivity_0p1pct.csv');

fprintf('\nMedian sensitivities (DGM)\n');
fprintf('T   = %.3f\n', median(SensTable.Sens_T_DGM, 'omitnan'));
fprintf('S   = %.3f\n', median(SensTable.Sens_S_DGM, 'omitnan'));
fprintf('Ice = %.3f\n', median(SensTable.Sens_Ice_DGM, 'omitnan'));
fprintf('W   = %.3f\n', median(SensTable.Sens_W_DGM, 'omitnan'));


function out = run_CAO_case(s_p, Tw, ice, u_10, lat, GEM_obs, DGM_obs)

    s_bg = 34.8;

    t_68 = (Tw+273.15-0.00025*273.15)/(1-0.00025)-273.15;
    p_sw_A = 0.824493 - 4.0899e-3*t_68 + 7.6438e-5*t_68^2 -8.2467e-7*t_68^3 + 5.3875e-9*t_68^4;
    p_sw_B = -5.72466e-3 + 1.0227e-4*t_68 - 1.6546e-6*t_68^2;
    p_sw_C = 4.8314e-4;
    p_w = 999.842594 + 6.793952e-2*t_68 - 9.09529e-3*t_68^2 + 1.001685e-4*t_68^3 - 1.120083e-6*t_68^4 + 6.536336e-9*t_68^5;
    p_sw = (p_w + p_sw_A*s_p + p_sw_B*s_p^(3/2) + p_sw_C*s_p)/1000;

    f_ow = s_p / s_bg;
    f_meteo = 1-f_ow;

    p_uhw = 1024.8;

    if ice >= 0.2
        MLD = 3.5337*s_p-91.066;
    else
        MLD = 3.6 * (p_uhw - p_sw*1000)^(-0.45)*(u_10)^(0.55);
    end

   A = 6.891e11; %m2, Area of the mixed layer in the CAO, rstudio
V = A*MLD; %m3, water volume

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
v_dep = 1.5; %cm s-1, Orbist et al., 2017 % 보완?
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
fr_THg = 4.0*1e-12*200.59*1e3; %g/m3, THg in the Makarov Basin, Eom et al., 2025       5.46
fr_ratio = 0.04834; % median of the CS, MB, and ESS;
fr_Hg2 = fr_THg*(1-fr_ratio);
fr_Hg0 = fr_THg*fr_ratio;
fr(1) = fr_vol*fr_Hg0;
fr(2) = fr_vol*fr_Hg2;

% Air-sea exchange                       
H = exp(-2404.3/(Tw+273.15)+6.92); % Dimensionless Henry's law constant
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
k_w = 0.25*(u_10)^2*(Sc_Hg/660)^(-0.5)*0.01/3600; %m s-1, mass transfer coefficient

%Particle - Kim 방법
Kp = 10^(6.74); % Cui et al., 2021, Overall Arctic Ocean
C_spm = 21.4*10^-9; % Phoebe J. Lam, Size-fractionated particle composition, Arctic Ocean
F_Hg2 = 1-(1/(1+Kp*C_spm)); % Particle fraction
V_set = 24.9*10^(-3)/24/3600*A; % g s-1, particle settling velocity, fall averaged without outliers

%Redox 
UV_B_in = 0.00015787*lat^2-0.026359*lat +1.1026; % Latitude 경험식 65 - 85, 140 - 205
UV_A_in = 0.01044*lat^2-1.8153*lat +79.325; % Latitude 경험식 65 - 85, 140 - 205
PAR_in = (0.0217*lat^2-3.8897*lat+178.7374)*10^6/24/3600/4.6; % Latitude 경험식 65 - 85, 140 - 205
I_trans = 0.028; % Median light transmittance, Katlein et al., 2019 (0.028 median) (0.011 median september)
UV_B = UV_B_in*(1-ice) + UV_B_in*ice*I_trans; %W m-2, intensity
UV_A = UV_A_in*(1-ice) + UV_A_in*ice*I_trans; %W m-2, intensity 
PAR = PAR_in*(1-ice) + PAR_in*ice*I_trans;  %W m-2, intensity

Chl = 0.0973*s_p-1.9444;

a_ph_443 = 10^(0.5877*log10(Chl)-1.5160); % Lewis and Arrigo,
DOC = -3.85 * s_p + 190; % Northern Chukchi Sea, Jung et al., 2022
a_g_443 = (DOC-55)/357; % CDOM absorption, Matsuoka et al., 2014
a_d_443 = 0.016; % Lewis and Arrigo,


a_tw_443 = 0.0145; %450 nm, Soerensen et al., 2010

a_443 = a_ph_443 + a_g_443 + a_d_443 + a_tw_443; 
Kd_490 = 0.0166+0.0733*Chl^0.6715; %(70 - 80oN)

Kd_380 = 3.313*(a_443)^2+0.887*a_443+0.03; %diffuse attenuation coefficient at 380 nm
Kd_340 = -1.965*(Kd_380)^2+1.979*Kd_380-0.009; %diffuse attenuation coefficient at 340 nm
Kd_305 = -7.07*(Kd_380)^2+4.424*Kd_380-0.008; %diffuse attenuation coefficient at 305 nm
Kd_PAR = 0.00665+0.874*Kd_490-0.00121*Kd_490^(-1); %diffuse attenuation coefficient at PAR

z_eu = log(100)/Kd_PAR;
D_irr = 16;
PAR_unit = PAR * 0.395;
Pb_opt = 2.0 * 10^(-8.7*10^(-3)*(MLD/2));% Lee et al., 2015, outperformed method -> Model 14 use equation (3) in Huot et al., 2013 if Chl-a > 0.1

NPP = 0.66125 * z_eu * Chl * D_irr * (PAR_unit/(PAR_unit+4.1))*Pb_opt;
pe_ratio = -0.0081 * Tw + 0.0668 * log(Chl/z_eu)+0.426;
OCRR_unit = NPP/z_eu * (1-pe_ratio); % mg m-3 day-12
OCRR = OCRR_unit / 1000  / 24 /12 ; % mol m-3 h-1



s0 = 29.0;
lambda_rhg = 0.05;

f_sal_rhg = 1 + lambda_rhg*(s_p - s0);
f_sal_rhg = max(0.90, min(1.10, f_sal_rhg));

RHg_UVA = ((-0.0274*(s_p)^2 + 1.6344*s_p - 23.752) / (1+Kp*C_spm)) * f_sal_rhg;


%RHg_UVA = (-0.0274*(s_p)^2+1.6344*(s_p)-23.752)/(1+Kp*C_spm); %Reducible Hg(II) fraction for UVA;                      %% 변경
RHg_UVB = 0.869/(1+Kp*C_spm);

UV_B_exp = 1.4; % W m-2, UV-B intensity in this study
UV_A_exp = 6.2; % W m-2, UV-A intensity in this study

kr_UVA = 0.084/UV_A_exp; %h-1 W-1 m2, light-normalized kr_UVA              %% 변경

kr_UVA_d = kr_UVA*UV_A*(1-exp(-Kd_340*MLD))/Kd_340/MLD; % depth-averaged kr_UVA
k_UVA(2) = kr_UVA_d*RHg_UVA; %Gross Hg(0) production rate (h-1) for UVA

kr_UVB = (0.8495*(RHg_UVA)^2-0.4843*RHg_UVA+0.2377)/UV_B_exp;

kr_UVB_d = kr_UVB*UV_B*(1-exp(-Kd_305*MLD))/Kd_305/MLD; % depth-averaged kr_UVB
k_UVB(2) = kr_UVB_d*RHg_UVB; %Gross Hg(0) production rate (h-1) for UVB

if MLD >=z_eu, kr_Bio_d = 86 * OCRR;
else kr_Bio_d = 86 * OCRR * MLD / z_eu;
end

k_Bio(2) = kr_Bio_d;

k_2 = k_UVB(2)+k_UVA(2)+k_Bio(2)*2;%+ k_Dark(2);%+k_PAR(2);

k(2) = k_2/3600;

UV_B_exp2 = 1.35;
UV_A_exp2 = 5.33;

s0 = 29.0;              % 기준 salinity, 데이터 중앙 근처
lambda_ko = 0.25;       % 시작값: 0.02~0.08 범위 탐색 추천

f_sal_ko = 1 + lambda_ko*(s0 - s_p);
f_sal_ko = max(0.85, min(1.15, f_sal_ko));   % 과도보정 방지

ko_UVA = ((0.0373*(s_p)^2 - 2.1285*s_p + 30.642) / UV_A_exp2) * f_sal_ko;


%ko_UVA = (0.0373*(s_p)^2-2.1285*s_p+30.642)/UV_A_exp2; % strong correlationm outlierX

%ko_UVA = 0.34/UV_A_exp2;
ko_UVB = 1.77/UV_B_exp2; %h-1 W-1 m2, light-normalized kr_UVB               %% 변경

ko_UVB_d = ko_UVB*UV_B*(1-exp(-Kd_305*MLD))/Kd_305/MLD; % depth-averaged kr_UVB
ko_UVA_d = ko_UVA*UV_A*(1-exp(-Kd_340*MLD))/Kd_340/MLD; % depth-averaged kr_UVA
k_UVB(1) = ko_UVB_d; %Gross Hg(0) production rate (h-1) for UVB
k_UVA(1) = ko_UVA_d; %Gross Hg(0) production rate (h-1) for UVA

k_Dark(1) = 1*10^(-7)*3600; % Soerensen et al., 2010; 

if MLD >=z_eu, k_Bio(1) = 140 * OCRR;
else k_Bio(1) = 140 * OCRR * MLD / z_eu;
end

k_1 = k_UVB(1)+k_UVA(1)+k_Dark(1)+k_Bio(1); %sum
k(1) = k_1/3600;

% Water column diffusion - Soerensen 방법 
Diff_d = 50/2;
g = 9.8;
r = 0.2; % empirical coefficient, Randelhoff et al., 
p_air = 1.3; % Arctic air density

Cd_oa = 1.5*10^(-3); %ocean - air drag coefficient, Andreas and Murphy, 1986
Cd_oi = 5*10^(-3); %ocean - ice drag coefficient, Fer et al., 2022 
tou_oa = (1-ice)*p_air*Cd_oa*(u_10)^2; % Calculation based on Loose et al., 2014
tou_oi = ice*p_sw*1000*Cd_oi*(0.015*u_10)^2; % Calculation based on Loose et al., 2014, 0.015 = Brenner et al., 2021
kappa = 0.41;

u_fric = sqrt((tou_oa+tou_oi)/(p_sw*1000));

e_diss = u_fric^3/(kappa*MLD);

N_bv =  g / (p_sw*1000) * (p_uhw - p_sw*1000) / Diff_d; % s-2, Brunt-Vaisala frequency
Dw = r * e_diss * (N_bv)^(-1);  %eddy diffusion coefficient

Cw(1) = 298 *10^(-15)*200.59*1000; % Subsurface Hg(0)
Cw(2) = 826 *10^(-15)*200.59*1000;% Subsurface Hg(II)

c_ss = solve_CAO_ss_direct( ...
    A,V,M_vol,MB,CS,ESS,CAO_nt,CAO_et, ...
    V_set,F_Hg2,Kp,k,H,k_w,ice,Cw,Dw,Diff_d,V_A, ...
    L_ess,L_cs_v,L_w_v,L_n,L_e,L_v_v,V_dep,air_ratio,P,fr);

c1_ss = c_ss(1)*1e9/200.59;
c2_ss = c_ss(2)*1e9/200.59;
c3_ss = c_ss(3)*1e9;

S0_Hg0toHgII = k(1)*V*c_ss(1);
S0_HgIItoHg0 = k(2)*V*c_ss(2);
AS = (1-ice)*k_w*(c_ss(1)-c_ss(3)/H)*A;
Net = S0_HgIItoHg0 - S0_Hg0toHgII;

out.Hg0_water = c1_ss;
out.HgII_water = c2_ss;
out.Hg0_air = c3_ss;
out.AS = AS;
out.Net = Net;
out.MLD = MLD;
out.Dw = Dw;

out.Red_UVB = k_UVB(2);
out.Red_UVA = k_UVA(2);
out.Red_Bio = k_Bio(2);
out.Ox_UVB = k_UVB(1);
out.Ox_UVA = k_UVA(1);
out.Ox_Dark = k_Dark(1);
out.Ox_Bio = k_Bio(1);
out.a_443 = a_443;
end