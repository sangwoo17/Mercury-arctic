clc; clear;

%% =========================
% 0) Load data and common subset
% =========================
data = readtable('inputs/DGM_Eval_CAO.csv');

ice_thr = 0.15;
s_min   = 26.46;

isIce0  = data.ICE > ice_thr;
isOpen0 = data.ICE <= ice_thr;

validIce  = isIce0  & data.Salinity >= (s_min + 1.1);
validOpen = isOpen0 & data.Salinity >= (s_min + 1.7);
useIdx = validIce | validOpen;

data_sub = data(useIdx,:);

Regime0 = strings(height(data_sub),1);
Regime0(data_sub.ICE > ice_thr)  = "ice";
Regime0(data_sub.ICE <= ice_thr) = "open";
data_sub.Regime0 = Regime0;

%% =========================
% 1) 5x5 phase-space settings
% =========================
sal_list = [0, -0.25, -0.50, -0.75, -1.00];
ice_mult_list = [1.00, 0.75, 0.50, 0.25, 0.00];

AllResults = table();
SummaryTable = table();

for ii = 1:numel(ice_mult_list)
    ice_mult = ice_mult_list(ii);

    for jj = 1:numel(sal_list)
        sal_anom = sal_list(jj);

        fprintf('Running phase scenario: ice x %.2f, sal %.2f psu\n', ice_mult, sal_anom);

        T = run_phase_scenario(data_sub, ice_thr, ice_mult, sal_anom);
        AllResults = [AllResults; T];

        S = summarize_phase(T, ice_mult, sal_anom);
        SummaryTable = [SummaryTable; S];
    end
end

writetable(AllResults, 'outputs/Hg_phase_space_all_5x5.csv');
writetable(SummaryTable, 'outputs/Hg_phase_space_summary_5x5.csv');

%% =========================
% 2) Baseline delta summary
% =========================
% baseline = ice=1.0, sal=0
baseIdx = abs(SummaryTable.IceMultiplier - 1.0) < 1e-12 & ...
          abs(SummaryTable.SalinityAnomaly - 0.0) < 1e-12;

Base = SummaryTable(baseIdx,:);

SummaryTable.d_Hg0  = SummaryTable.Hg0_water_median - Base.Hg0_water_median;
SummaryTable.d_HgII = SummaryTable.HgII_water_median - Base.HgII_water_median;
SummaryTable.d_Net  = SummaryTable.Net_median - Base.Net_median;
SummaryTable.d_AS   = SummaryTable.AS_median - Base.AS_median;

SummaryTable.pct_Hg0  = 100 * SummaryTable.d_Hg0  / Base.Hg0_water_median;
SummaryTable.pct_HgII = 100 * SummaryTable.d_HgII / Base.HgII_water_median;
SummaryTable.pct_Net  = 100 * SummaryTable.d_Net  / Base.Net_median;
SummaryTable.pct_AS   = 100 * SummaryTable.d_AS   / Base.AS_median;

% dominance metric
SummaryTable.DomScore = abs(SummaryTable.d_Net) - abs(SummaryTable.d_AS);
% DomScore > 0 : chemistry-dominated
% DomScore < 0 : exchange-dominated

writetable(SummaryTable, 'outputs/Hg_phase_space_summary_5x5_delta.csv');

disp('Saved:')
disp(' - outputs/Hg_phase_space_all_5x5.csv')
disp(' - outputs/Hg_phase_space_summary_5x5.csv')
disp(' - outputs/Hg_phase_space_summary_5x5_delta.csv')

%% =========================
% Local function: phase scenario
% =========================
function ResultTable = run_phase_scenario(data_sub, ice_thr, ice_mult, sal_anom)

ResultTable = table();

for i = 1:height(data_sub)

    s_bg   = 34.8;
    s_p0   = data_sub.Salinity(i);
    Tw     = data_sub.Temperature(i);
    ice0   = data_sub.ICE(i);
    u_10   = data_sub.WS(i);
    lat    = data_sub.Lat(i);
    GEM_obs = data_sub.GEM(i);
    DGM_obs = data_sub.DGM(i);
    regime0 = data_sub.Regime0(i);

    % Phase-space perturbation
    ice = ice0 * ice_mult;
    s_p = s_p0 + sal_anom;

    if s_p < 26.46
        continue
    end

    if ice > ice_thr
        regime_new = "ice";
    else
        regime_new = "open";
    end

    % =========================
    % Model body
    % =========================
    t_68 = (Tw+273.15-0.00025*273.15)/(1-0.00025)-273.15;
    p_sw_A = 0.824493 - 4.0899e-3*t_68 + 7.6438e-5*t_68^2 - 8.2467e-7*t_68^3 + 5.3875e-9*t_68^4;
    p_sw_B = -5.72466e-3 + 1.0227e-4*t_68 - 1.6546e-6*t_68^2;
    p_sw_C = 4.8314e-4;
    p_w = 999.842594 + 6.793952e-2*t_68 - 9.09529e-3*t_68^2 + 1.001685e-4*t_68^3 ...
        - 1.120083e-6*t_68^4 + 6.536336e-9*t_68^5;
    p_sw = (p_w + p_sw_A*s_p + p_sw_B*s_p^(3/2) + p_sw_C*s_p)/1000;

    f_ow = s_p/s_bg;
    f_meteo = 1-f_ow;
    p_uhw = 1024.8;

    if ice > ice_thr
        MLD = 3.5337*s_p - 91.066;
    else
        MLD = 3.6*(p_uhw - p_sw*1000)^(-0.45)*(u_10)^(0.55);
    end

    A = 6.891e11;
    V = A*MLD;
    ABL = 236;
    V_A = A*ABL;
    air_ratio = 0.00733;

    % Atmospheric transport
    L_ess_wind = 3.1475; L_ess_dist = 1038.32e3;
    L_ess_v = L_ess_wind * L_ess_dist * ABL;
    air_Hg_ess = 1.05; L_ess(1) = L_ess_v*air_Hg_ess*1e-9;

    L_cs_wind = -2.5303; L_cs_dist = 610.01e3;
    L_cs_v = L_cs_wind * L_cs_dist * ABL;

    L_w_wind = -4.5416; L_w_dist = 500e3;
    L_w_v = L_w_wind * L_w_dist * ABL;

    L_n_wind = -2.4825; L_n_dist = 1105e3;
    L_n_v = L_n_wind * L_n_dist * ABL;
    air_Hg_n = 1.176829; L_n(1) = L_n_v*air_Hg_n*1e-9;

    L_e_wind = -2.5448; L_e_dist = 500e3;
    L_e_v = L_e_wind * L_e_dist * ABL;
    air_Hg_e = 1.144374; L_e(1) = L_e_v*air_Hg_e*1e-9;

    L_v_v = L_ess_v + L_cs_v + L_w_v - L_e_v - L_n_v;

    % Deposition
    v_dep = 1.5;
    V_dep = v_dep*0.01*A;
    precip = 0.001061/24/3600*A;
    precip_thg = 3.1e-12*200.59*1000;
    P = (1-ice)*precip*precip_thg;

    % Meltwater
    ice_b = 176.07e9; ice_a = 278.32e9;
    M_vol = (ice_b-ice_a)/30/86400;

    % Ocean current
    MB_v = 0.022997; MB_d = 500377.2;
    MB_vol = MB_v*MLD*MB_d;
    MB_THg = 1.6e-12*200.59*1e3; MB_ratio = 0.07535;
    MB_Hg2 = MB_THg*(1-MB_ratio); MB_Hg0 = MB_THg*MB_ratio;
    MB(1) = MB_vol*MB_Hg0*(1-f_meteo); MB(2) = MB_vol*MB_Hg2*(1-f_meteo);

    CS_v = 0.035393; CS_d = 610011.6;
    CS_vol = CS_v*MLD*CS_d;
    CS_THg = 1.7e-12*200.59*1e3; CS_ratio = 0.02647;
    CS_Hg2 = CS_THg*(1-CS_ratio); CS_Hg0 = CS_THg*CS_ratio;
    CS(1) = CS_vol*CS_Hg0*(1-f_meteo); CS(2) = CS_vol*CS_Hg2*(1-f_meteo);

    ESS_v = 0.011518; ESS_d = 1038317.6;
    ESS_vol = ESS_v*MLD*ESS_d;
    ESS_THg = 1.0e-12*200.59*1e3; ESS_ratio = 0.04834;
    ESS_Hg2 = ESS_THg*(1-ESS_ratio); ESS_Hg0 = ESS_THg*ESS_ratio;
    ESS(1) = ESS_vol*ESS_Hg0*(1-f_meteo); ESS(2) = ESS_vol*ESS_Hg2*(1-f_meteo);

    CAO_nv = 0.016198; CAO_nd = 1104564.6;
    CAO_nt = CAO_nv*MLD*CAO_nd;

    CAO_ev = 0.029777; CAO_ed = 500377.2;
    CAO_et = CAO_ev*MLD*CAO_ed;

    % Freshwater
    fr_vol = (ESS_vol+CS_vol+MB_vol)*f_meteo;
    fr_THg = 4.0e-12*200.59*1e3; fr_ratio = 0.04834;
    fr_Hg2 = fr_THg*(1-fr_ratio); fr_Hg0 = fr_THg*fr_ratio;
    fr(1) = fr_vol*fr_Hg0; fr(2) = fr_vol*fr_Hg2;

    % Air-sea exchange
    H = exp(-2404.3/(Tw+273.15)+6.92);
    n_fw = 4.2844e-5 + (0.157*(Tw+64.993)^2 - 91.296)^(-1);
    n_sw_A = 1.541 + 1.998e-2*Tw - 9.52e-5*Tw^2;
    n_sw_B = 7.974 - 7.561e-2*Tw + 4.724e-4*Tw^2;
    n_sw = 10*n_fw*(1+n_sw_A*(s_p)/1000+n_sw_B*(s_p/1000)^2);
    pi_w = 2.26; M_w = 18.01; V_b = 12.74;
    D_sw_Hg = 7.4e-10*(273.15+Tw)*sqrt(pi_w*M_w)/n_sw/V_b^0.6;
    v_sw = n_sw/p_sw;
    Sc_Hg = v_sw/D_sw_Hg;
    k_w = 0.25*(u_10)^2*(Sc_Hg/660)^(-0.5)*0.01/3600;

    % Particle
    Kp = 10^(6.74); C_spm = 21.4e-9;
    F_Hg2 = 1-(1/(1+Kp*C_spm));
    V_set = 24.9e-3/24/3600*A;

    % Redox
    UV_B_in = 0.00015787*lat^2 - 0.026359*lat + 1.1026;
    UV_A_in = 0.01044*lat^2 - 1.8153*lat + 79.325;
    PAR_in  = (0.0217*lat^2 - 3.8897*lat + 178.7374)*1e6/24/3600/4.6;
    I_trans = 0.028;
    UV_B = UV_B_in*(1-ice) + UV_B_in*ice*I_trans;
    UV_A = UV_A_in*(1-ice) + UV_A_in*ice*I_trans;
    PAR  = PAR_in*(1-ice) + PAR_in*ice*I_trans;

    Chl = 0.0973*s_p - 1.9444;
    a_ph_443 = 10^(0.5877*log10(Chl)-1.5160);
    DOC = -3.85*s_p + 190;
    a_g_443 = (DOC-55)/357;
    a_d_443 = 0.016; a_tw_443 = 0.0145;
    a_443 = a_ph_443 + a_g_443 + a_d_443 + a_tw_443;

    Kd_490 = 0.0166 + 0.0733*Chl^0.6715;
    Kd_380 = 3.313*(a_443)^2 + 0.887*a_443 + 0.03;
    Kd_340 = -1.965*(Kd_380)^2 + 1.979*Kd_380 - 0.009;
    Kd_305 = -7.07*(Kd_380)^2 + 4.424*Kd_380 - 0.008;
    Kd_PAR = 0.00665 + 0.874*Kd_490 - 0.00121*Kd_490^(-1);

    z_eu = log(100)/Kd_PAR;
    D_irr = 16;
    PAR_unit = PAR*0.395;
    Pb_opt = 2.0 * 10^(-8.7e-3*(MLD/2));

    NPP = 0.66125*z_eu*Chl*D_irr*(PAR_unit/(PAR_unit+4.1))*Pb_opt;
    pe_ratio = -0.0081*Tw + 0.0668*log(Chl/z_eu) + 0.426;
    OCRR_unit = NPP/z_eu*(1-pe_ratio);
    OCRR = OCRR_unit/1000/24/12;

    s0 = 29.0;
    lambda_rhg = 0.05;
    f_sal_rhg = 1 + lambda_rhg*(s_p - s0);
    f_sal_rhg = max(0.90, min(1.10, f_sal_rhg));
    RHg_UVA = ((-0.0274*s_p^2 + 1.6344*s_p - 23.752)/(1+Kp*C_spm))*f_sal_rhg;
    RHg_UVB = 0.869/(1+Kp*C_spm);

    UV_B_exp = 1.4; UV_A_exp = 6.2;
    kr_UVA = 0.084/UV_A_exp;
    kr_UVA_d = kr_UVA*UV_A*(1-exp(-Kd_340*MLD))/Kd_340/MLD;
    k_UVA(2) = kr_UVA_d*RHg_UVA;

    kr_UVB = (0.8495*(RHg_UVA)^2 - 0.4843*RHg_UVA + 0.2377)/UV_B_exp;
    kr_UVB_d = kr_UVB*UV_B*(1-exp(-Kd_305*MLD))/Kd_305/MLD;
    k_UVB(2) = kr_UVB_d*RHg_UVB;

    if MLD >= z_eu
        kr_Bio_d = 86*OCRR;
    else
        kr_Bio_d = 86*OCRR*MLD/z_eu;
    end
    k_Bio(2) = kr_Bio_d;
    k_2 = k_UVB(2)+k_UVA(2)+k_Bio(2)*2;
    k(2) = k_2/3600;

    UV_B_exp2 = 1.35; UV_A_exp2 = 5.33;
    lambda_ko = 0.25;
    f_sal_ko = 1 + lambda_ko*(s0 - s_p);
    f_sal_ko = max(0.85, min(1.15, f_sal_ko));
    ko_UVA = ((0.0373*s_p^2 - 2.1285*s_p + 30.642)/UV_A_exp2)*f_sal_ko;
    ko_UVB = 1.77/UV_B_exp2;

    ko_UVB_d = ko_UVB*UV_B*(1-exp(-Kd_305*MLD))/Kd_305/MLD;
    ko_UVA_d = ko_UVA*UV_A*(1-exp(-Kd_340*MLD))/Kd_340/MLD;
    k_UVB(1) = ko_UVB_d;
    k_UVA(1) = ko_UVA_d;
    k_Dark(1) = 1e-7*3600;

    if MLD >= z_eu
        k_Bio(1) = 140*OCRR;
    else
        k_Bio(1) = 140*OCRR*MLD/z_eu;
    end
    k_1 = k_UVB(1)+k_UVA(1)+k_Dark(1)+k_Bio(1);
    k(1) = k_1/3600;

    % Diffusion
    Diff_d = 25;
    g = 9.8; r = 0.2; p_air = 1.3;
    Cd_oa = 1.5e-3; Cd_oi = 5e-3;
    tou_oa = (1-ice)*p_air*Cd_oa*(u_10)^2;
    tou_oi = ice*p_sw*1000*Cd_oi*(0.015*u_10)^2;
    kappa = 0.41;

    u_fric = sqrt((tou_oa+tou_oi)/(p_sw*1000));
    e_diss = u_fric^3/(kappa*MLD);
    N_bv = g/(p_sw*1000)*(p_uhw - p_sw*1000)/Diff_d;
    Dw = r*e_diss*(N_bv)^(-1);

    Cw(1) = 298e-15*200.59*1000;
    Cw(2) = 826e-15*200.59*1000;

    c_ss = solve_CAO_ss_direct( ...
        A,V,M_vol,MB,CS,ESS,CAO_nt,CAO_et, ...
        V_set,F_Hg2,Kp,k,H,k_w,ice,Cw,Dw,Diff_d,V_A, ...
        L_ess,L_cs_v,L_w_v,L_n,L_e,L_v_v,V_dep,air_ratio,P,fr);

    c1_ss = c_ss(1)*1e9/200.59;
    c2_ss = c_ss(2)*1e9/200.59;
    c3_ss = c_ss(3)*1e9;

    S0_HgIISink = V_set*c_ss(2)*(1-F_Hg2)*Kp*1e-6;
    S0_Hg0toHgII = k(1)*V*c_ss(1);
    S0_HgIItoHg0 = k(2)*V*c_ss(2);
    AS = (1-ice)*k_w*(c_ss(1)-c_ss(3)/H)*A;

    Diff_Hg0 = A*Dw*(Cw(1)-c_ss(1))/Diff_d;
    Diff_HgII = A*Dw*(Cw(2)-c_ss(2)*(1-F_Hg2))/Diff_d;
    Current_Hg0 = MB(1)+CS(1)+ESS(1);
    Current_HgII = MB(2)+CS(2)+ESS(2);
    Net = S0_HgIItoHg0 - S0_Hg0toHgII;

    f_UVB_2 = k_UVB(2); f_UVA_2 = k_UVA(2); f_Bio_2 = k_Bio(2);
    f_UVB_1 = k_UVB(1); f_UVA_1 = k_UVA(1); f_Dark_1 = k_Dark(1); f_Bio_1 = k_Bio(1);

    NewRow = table( ...
        i, string(regime0), string(regime_new), ice_mult, sal_anom, ...
        Tw, s_p0, s_p, ice0, ice, u_10, ...
        c1_ss, DGM_obs, c2_ss, c3_ss, GEM_obs, MLD, Dw, ...
        S0_HgIISink, AS, Net, Diff_Hg0, Current_Hg0, fr(1), ...
        Diff_HgII, Current_HgII, fr(2), ...
        f_UVB_2, f_UVA_2, f_Bio_2, f_UVB_1, f_UVA_1, f_Dark_1, f_Bio_1, ...
        'VariableNames', { ...
        'RowID','Regime0','RegimeNew','IceMultiplier','SalinityAnomaly', ...
        'Temperature','Salinity0','Salinity','SeaIce0','SeaIce','WindSpeed', ...
        'Hg0_water','DGM_obs','HgII_water','Hg0_air','GEM_obs','MLD','Dw', ...
        'P_sink','AS','Net','Diff_Hg0','Current_Hg0','fr_Hg0', ...
        'Diff_HgII','Current_HgII','fr_HgII', ...
        'Red_UVB','Red_UVA','Red_Bio','Ox_UVB','Ox_UVA','Ox_Dark','Ox_Bio'} ...
        );

    ResultTable = [ResultTable; NewRow];
end
end

%% =========================
% Local function: summary
% =========================
function S = summarize_phase(T, ice_mult, sal_anom)

vars = {'Hg0_water','HgII_water','MLD','AS','Net','Red_UVA','Ox_UVA','SeaIce','Salinity'};

S = table(ice_mult, sal_anom, height(T), ...
    'VariableNames', {'IceMultiplier','SalinityAnomaly','N'});

for k = 1:numel(vars)
    v = vars{k};
    S.([v '_median']) = median(T.(v), 'omitnan');
    S.([v '_mean'])   = mean(T.(v), 'omitnan');
end

S.n_shift = sum(T.Regime0 ~= T.RegimeNew);
end

%% =========================
% Local function: steady-state solver
% =========================
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

M(1,1) = -k(1) + M_vol/V - CAO_nt/V - CAO_et/V - kwA/V - diff0/V;
M(1,2) =  k(2);
M(1,3) =  kwA/(H*V);
rhs(1) = -InHg0/V;

M(2,1) =  k(1);
M(2,2) = -k(2) - sinkHgII/V - CAO_nt/V - CAO_et/V - diff2_out/V;
M(2,3) =  (1-ice)*depCoeff/V;
rhs(2) = -InHgII/V;

M(3,1) =  kwA/V_A;
M(3,2) =  0;
M(3,3) = (L_w_v + L_cs_v - L_v_v)/V_A - kwA/(H*V_A) - depCoeff/V_A;
rhs(3) = -InAir/V_A;

c_ss = M \ rhs;
end



