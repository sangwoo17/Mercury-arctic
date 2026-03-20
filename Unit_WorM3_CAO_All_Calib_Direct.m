clc; clear; 

data = readtable('DGM_Valid_CAO.csv');
ResultTable = table();

for i = 1:height(data)

% Boundary conditions

s_bg = 34.8;
s_p = data.Salinity(i); % or s_p_raw

Tw= data.Temperature(i); % degree C, in situ seawater temperature, ITS-90 기준   

ice = data.ICE(i);                           
u_10 = data.WS(i); % wind speed
lat = data.Lat(i); % latitude

GEM_obs = data.GEM(i);
DGM_obs = data.DGM(i);

t_68 = (Tw+273.15-0.00025*273.15)/(1-0.00025)-273.15; % degree C, in situ seawater temperature, IPTS-68 기준
p_sw_A = 0.824493 - 4.0899*10^(-3) * t_68 + 7.6438 * 10^(-5) * (t_68)^2 -8.2467 * 10^(-7) *(t_68)^3 + 5.3875*10^(-9) * (t_68)^4;
p_sw_B = -5.72466*10^(-3)+1.0227*10^(-4)*t_68-1.6546*10^(-6)*(t_68)^2;
p_sw_C = 4.8314*10^(-4);
p_w = 999.842594 + 6.793952*10^(-2)*t_68 - 9.09529*10^(-3)*(t_68)^2 + 1.001685*10^(-4)*(t_68)^3 - 1.120083*10^(-6)*(t_68)^4 + 6.536336*10^(-9)*(t_68)^5;
p_sw = (p_w +  p_sw_A * s_p + p_sw_B * s_p^(3/2) + p_sw_C * s_p)/1000;
f_ow = s_p / s_bg;
f_meteo = 1-f_ow;

p_wml = 1024.6;



if ice >= 0.15
    % Ice-covered ocean
    MLD = 3.5337*s_p-91.066;

else
    % Open water
    MLD = 3.6 * (p_wml - p_sw*1000)^(-0.45)*(u_10)^(0.55); %m, mixed layer depth of the Makarov Basin
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
P = precip * precip_thg; % g s-1, wet deposition

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
Kp = 10^(6.83); % Cui et al., 2021, Overall Arctic Ocean
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

RHg_UVA = ((-0.0274*(s_p)^2 + 1.6344*s_p - 23.752) / (1+Kp*C_spm));
RHg_UVB = 0.869/(1+Kp*C_spm);

UV_B_exp = 1.4; % W m-2, UV-B intensity in this study
UV_A_exp = 6.2; % W m-2, UV-A intensity in this study

kr_UVA = (0.0088*s_p^2-0.513*s_p+7.5412)/UV_A_exp; %h-1 W-1 m2, light-normalized kr_UVA              %% 변경
kr_UVA_d = kr_UVA*UV_A*(1-exp(-Kd_340*MLD))/Kd_340/MLD; % depth-averaged kr_UVA
k_UVA(2) = kr_UVA_d*RHg_UVA; %Gross Hg(0) production rate (h-1) for UVA

kr_UVB = 0.207/UV_B_exp; % CAO only
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

ko_UVA = (0.340)/UV_A_exp2;
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

N_bv =  g / (p_sw*1000) * (p_wml - p_sw*1000) / Diff_d; % s-2, Brunt-Vaisala frequency
Dw = r * e_diss * (N_bv)^(-1) ;  %eddy diffusion coefficient

Cw(1) = 298 *10^(-15)*200.59*1000; % Subsurface Hg(0)
Cw(2) = 826 *10^(-15)*200.59*1000;% Subsurface Hg(II)

c_ss = solve_CAO_ss_direct( ...
    A,V,M_vol,MB,CS,ESS,CAO_nt,CAO_et, ...
    V_set,F_Hg2,Kp,k,H,k_w,ice,Cw,Dw,Diff_d,V_A, ...
    L_ess,L_cs_v,L_w_v,L_n,L_e,L_v_v,V_dep,air_ratio,P,fr);

c1_ss = c_ss(1)*10^9/200.59;
c2_ss = c_ss(2)*10^9/200.59;
c3_ss = c_ss(3)*10^9;

S0_HgIISink = V_set*c_ss(2)*(1-F_Hg2)*Kp*10^(-6);
S0_Hg0toHgII = k(1)*V*c_ss(1);
S0_HgIItoHg0 = k(2)*V*c_ss(2);
AS = (1-ice)*k_w*(c_ss(1)-c_ss(3)/H)*A;
fr_HgII = fr(2);
fr_Hg0 = fr(1);
Diff_Hg0 = A*Dw*(Cw(1)-c_ss(1))/Diff_d;
Diff_HgII = A*Dw*(Cw(2)-c_ss(2)*(1-F_Hg2))/Diff_d;
Current_Hg0 = MB(1)+CS(1)+ESS(1);
Current_HgII = MB(2)+CS(2)+ESS(2);


NSE_num = (DGM_obs - c1_ss)^2;
NSE_den = (DGM_obs - 0.12056)^2;
PBIAS_num = (DGM_obs - c1_ss)*100;
PBIAS_den = DGM_obs;
Net = S0_HgIItoHg0 - S0_Hg0toHgII;
f_UVB_2 = k_UVB(2);
f_UVA_2 = k_UVA(2);
f_Bio_2 = k_Bio(2);
f_UVB_1 = k_UVB(1);
f_UVA_1 = k_UVA(1);
f_Dark_1 = k_Dark(1);
f_Bio_1 = k_Bio(1);


    NewRow = table( ...
        i, Tw, s_p, ice, u_10, ...
        c1_ss, DGM_obs, c2_ss, c3_ss, GEM_obs, MLD, Dw, ...
        S0_HgIISink,AS ,Net, Diff_Hg0,  Current_Hg0,  fr_Hg0, ...
         Diff_HgII,Current_HgII, fr_HgII, f_UVB_2, f_UVA_2, f_Bio_2, f_UVB_1, f_UVA_1, f_Dark_1, f_Bio_1, a_443, ...
        NSE_num, NSE_den, PBIAS_num, PBIAS_den, ...
        'VariableNames', { ...
        'RowID', 'Temperature', 'Salinity', 'SeaIce', 'WindSpeed', ...
        'Hg0_water', 'DGM_obs', 'HgII_water', 'Hg0_air', 'GEM_obs', ...
        'MLD', 'Dw', 'P_sink', 'AS', 'Net', 'Diff_Hg0',  'Current_Hg0', 'fr_Hg0', ...
        'Diff_HgII','Current_HgII', 'fr_HgII', 'Red_UVB', 'Red_UVA', 'Red_Bio', 'Ox_UVB', 'Ox_UVA', 'Ox_Dark', 'Ox_Bio', 'a_443', 'NSE_num', 'NSE_den', 'PBIAS_num', 'PBIAS_den'} ...
    );

    ResultTable = [ResultTable; NewRow];

end

writetable(ResultTable, 'Hg_budget_summary.csv');



%% =========================
% Model performance metrics
% No toolbox required
% =========================
obs = ResultTable.DGM_obs;
sim = ResultTable.Hg0_water;

% NaN 제거
valid = isfinite(obs) & isfinite(sim);
obs = obs(valid);
sim = sim(valid);

n = numel(obs);

% 1) Linear regression: sim = slope * obs + intercept
p = polyfit(obs, sim, 1);
slope = p(1);
intercept = p(2);

% 회귀선 예측값
sim_fit = polyval(p, obs);

% 2) R^2 (regression-based)
SS_res_fit = sum((sim - sim_fit).^2);
SS_tot_fit = sum((sim - mean(sim)).^2);
R2_reg = 1 - SS_res_fit / SS_tot_fit;

% 3) Pearson r and r^2 (manual, no toolbox)
obs_mean = mean(obs);
sim_mean = mean(sim);

num = sum((obs - obs_mean) .* (sim - sim_mean));
den = sqrt(sum((obs - obs_mean).^2) * sum((sim - sim_mean).^2));

r = num / den;
R2_corr = r^2;

% 4) RMSE
RMSE = sqrt(mean((sim - obs).^2));

% 5) PBIAS
PBIAS = 100 * sum(obs - sim) / sum(obs);

% 6) Mean Bias
Bias = mean(sim - obs);

% 출력
fprintf('\n===== Model Performance Metrics =====\n');
fprintf('N             = %d\n', n);
fprintf('R^2 (reg)     = %.4f\n', R2_reg);
fprintf('R^2 (corr^2)  = %.4f\n', R2_corr);
fprintf('r             = %.4f\n', r);
fprintf('RMSE          = %.4f\n', RMSE);
fprintf('Slope         = %.4f\n', slope);
fprintf('Y-intercept   = %.4f\n', intercept);
fprintf('PBIAS (%%)     = %.2f\n', PBIAS);
fprintf('Mean Bias     = %.4f\n', Bias);

% summary table 저장
MetricsTable = table(n, R2_reg, R2_corr, r, RMSE, slope, intercept, PBIAS, Bias, ...
    'VariableNames', {'N','R2_reg','R2_corr','r','RMSE','Slope','Intercept','PBIAS_percent','MeanBias'});

writetable(MetricsTable, 'Hg_model_metrics.csv');


function c_ss = solve_CAO_ss_direct(A,V,M_vol,MB,CS,ESS,CAO_nt,CAO_et,...
    V_set,F_Hg2,Kp,k,H,k_w,ice,Cw,Dw,Diff_d,V_A,...
    L_ess,L_cs_v,L_w_v,L_n,L_e,L_v_v,V_dep,air_ratio,P,fr)

kwA = (1-ice)*k_w*A;
diff0 = A*Dw/Diff_d;
diff2_in = A*Dw/Diff_d;
diff2_out = A*Dw*(1-F_Hg2)/Diff_d;
sinkHgII = V_set*(1-F_Hg2)*Kp*1e-6;
depCoeff = V_dep*air_ratio;

InHg0  = MB(1) + CS(1) + ESS(1) + fr(1) + diff0*Cw(1);
InHgII = MB(2) + CS(2) + ESS(2) + fr(2) + diff2_in*Cw(2) + (1-ice)*P;
InAir  = L_ess(1) - L_n(1) - L_e(1) - P;

M = zeros(3,3);
rhs = zeros(3,1);

% Hg(0) water
M(1,1) = -k(1) + M_vol/V - CAO_nt/V - CAO_et/V - kwA/V - diff0/V;
M(1,2) =  k(2);
M(1,3) =  kwA/(H*V);
rhs(1) = -InHg0/V;

% Hg(II) water
M(2,1) =  k(1);
M(2,2) = -k(2) - sinkHgII/V - CAO_nt/V - CAO_et/V - diff2_out/V;
M(2,3) =  (1-ice)*depCoeff/V;
rhs(2) = -InHgII/V;

% Hg(0) air
M(3,1) =  kwA/V_A;
M(3,2) =  0;
M(3,3) = (L_w_v + L_cs_v - L_v_v)/V_A - kwA/(H*V_A) - depCoeff/V_A;
rhs(3) = -InAir/V_A;

c_ss = M \ rhs;
end



