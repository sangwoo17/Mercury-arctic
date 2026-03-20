from __future__ import annotations

import numpy as np
import pandas as pd


def solve_cao_ss_direct(
    A,
    V,
    M_vol,
    MB,
    CS,
    ESS,
    cao_nt,
    cao_et,
    V_set,
    F_Hg2,
    Kp,
    k,
    H,
    k_w,
    ice,
    Cw,
    Dw,
    diff_d,
    V_A,
    L_ess,
    L_cs_v,
    L_w_v,
    L_n,
    L_e,
    L_v_v,
    V_dep,
    air_ratio,
    P,
    fr,
):
    kwA = (1 - ice) * k_w * A
    diff0 = A * Dw / diff_d
    diff2_in = A * Dw / diff_d
    diff2_out = A * Dw * (1 - F_Hg2) / diff_d
    sink_hgii = V_set * (1 - F_Hg2) * Kp * 1e-6
    dep_coeff = V_dep * air_ratio

    in_hg0 = MB[0] + CS[0] + ESS[0] + fr[0] + diff0 * Cw[0]
    in_hgii = MB[1] + CS[1] + ESS[1] + fr[1] + diff2_in * Cw[1] + (1 - ice) * P
    in_air = L_ess[0] - L_n[0] - L_e[0] - P

    M = np.zeros((3, 3), dtype=float)
    rhs = np.zeros(3, dtype=float)

    M[0, 0] = -k[0] + M_vol / V - cao_nt / V - cao_et / V - kwA / V - diff0 / V
    M[0, 1] = k[1]
    M[0, 2] = kwA / (H * V)
    rhs[0] = -in_hg0 / V

    M[1, 0] = k[0]
    M[1, 1] = -k[1] - sink_hgii / V - cao_nt / V - cao_et / V - diff2_out / V
    M[1, 2] = (1 - ice) * dep_coeff / V
    rhs[1] = -in_hgii / V

    M[2, 0] = kwA / V_A
    M[2, 1] = 0.0
    M[2, 2] = (L_w_v + L_cs_v - L_v_v) / V_A - kwA / (H * V_A) - dep_coeff / V_A
    rhs[2] = -in_air / V_A

    return np.linalg.solve(M, rhs)


def run_model(
    input_csv: str = "DGM_Valid_CAO.csv",
    output_summary_csv: str = "Hg_budget_summary.csv",
    output_metrics_csv: str = "Hg_model_metrics.csv",
):
    data = pd.read_csv(input_csv)
    rows = []

    for i, row in data.iterrows():
        s_bg = 34.8
        s_p = row["Salinity"]
        Tw = row["Temperature"]
        ice = row["ICE"]
        u_10 = row["WS"]
        lat = row["Lat"]

        GEM_obs = row["GEM"]
        DGM_obs = row["DGM"]

        t_68 = (Tw + 273.15 - 0.00025 * 273.15) / (1 - 0.00025) - 273.15
        p_sw_A = (
            0.824493
            - 4.0899e-3 * t_68
            + 7.6438e-5 * t_68**2
            - 8.2467e-7 * t_68**3
            + 5.3875e-9 * t_68**4
        )
        p_sw_B = -5.72466e-3 + 1.0227e-4 * t_68 - 1.6546e-6 * t_68**2
        p_sw_C = 4.8314e-4
        p_w = (
            999.842594
            + 6.793952e-2 * t_68
            - 9.09529e-3 * t_68**2
            + 1.001685e-4 * t_68**3
            - 1.120083e-6 * t_68**4
            + 6.536336e-9 * t_68**5
        )
        p_sw = (p_w + p_sw_A * s_p + p_sw_B * s_p ** (3 / 2) + p_sw_C * s_p) / 1000
        f_ow = s_p / s_bg
        f_meteo = 1 - f_ow

        p_wml = 1024.6

        if ice >= 0.15:
            MLD = 3.5337 * s_p - 91.066
        else:
            MLD = 3.6 * (p_wml - p_sw * 1000) ** (-0.45) * u_10**0.55

        A = 6.891e11
        V = A * MLD

        ABL = 236
        V_A = A * ABL
        air_ratio = 0.00733

        L_ess_wind = 3.1475
        L_ess_dist = 1038.32 * 1000
        L_ess_v = L_ess_wind * L_ess_dist * ABL
        air_Hg_ess = 1.05
        L_ess = np.array([L_ess_v * air_Hg_ess * 1e-9], dtype=float)

        L_cs_wind = -2.5303
        L_cs_dist = 610.01 * 1000
        L_cs_v = L_cs_wind * L_cs_dist * ABL

        L_w_wind = -4.5416
        L_w_dist = 500 * 1000
        L_w_v = L_w_wind * L_w_dist * ABL

        L_n_wind = -2.4825
        L_n_dist = 1105 * 1000
        L_n_v = L_n_wind * L_n_dist * ABL
        air_Hg_n = 1.176829
        L_n = np.array([L_n_v * air_Hg_n * 1e-9], dtype=float)

        L_e_wind = -2.5448
        L_e_dist = 500 * 1000
        L_e_v = L_e_wind * L_e_dist * ABL
        air_Hg_e = 1.144374
        L_e = np.array([L_e_v * air_Hg_e * 1e-9], dtype=float)

        L_v_v = L_ess_v + L_cs_v + L_w_v - L_e_v - L_n_v

        v_dep = 1.5
        V_dep = v_dep * 0.01 * A

        precip = 0.001061 / 24 / 3600 * A
        precip_thg = 3.1e-12 * 200.59 * 1000
        P = precip * precip_thg

        ice_b = 176.07e9
        ice_a = 278.32e9
        M_vol = (ice_b - ice_a) / 30 / 86400

        MB_v = 0.022997
        MB_d = 500377.2
        MB_vol = MB_v * MLD * MB_d
        MB_THg = 1.6e-12 * 200.59 * 1e3
        MB_ratio = 0.07535
        MB_Hg2 = MB_THg * (1 - MB_ratio)
        MB_Hg0 = MB_THg * MB_ratio
        MB = np.array([MB_vol * MB_Hg0 * (1 - f_meteo), MB_vol * MB_Hg2 * (1 - f_meteo)], dtype=float)

        CS_v = 0.035393
        CS_d = 610011.6
        CS_vol = CS_v * MLD * CS_d
        CS_THg = 1.7e-12 * 200.59 * 1e3
        CS_ratio = 0.02647
        CS_Hg2 = CS_THg * (1 - CS_ratio)
        CS_Hg0 = CS_THg * CS_ratio
        CS = np.array([CS_vol * CS_Hg0 * (1 - f_meteo), CS_vol * CS_Hg2 * (1 - f_meteo)], dtype=float)

        ESS_v = 0.011518
        ESS_d = 1038317.6
        ESS_vol = ESS_v * MLD * ESS_d
        ESS_THg = 1.0e-12 * 200.59 * 1e3
        ESS_ratio = 0.04834
        ESS_Hg2 = ESS_THg * (1 - ESS_ratio)
        ESS_Hg0 = ESS_THg * ESS_ratio
        ESS = np.array([ESS_vol * ESS_Hg0 * (1 - f_meteo), ESS_vol * ESS_Hg2 * (1 - f_meteo)], dtype=float)

        CAO_nv = 0.016198
        CAO_nd = 1104564.6
        CAO_nt = CAO_nv * MLD * CAO_nd

        CAO_ev = 0.029777
        CAO_ed = 500377.2
        CAO_et = CAO_ev * MLD * CAO_ed

        fr_vol = (ESS_vol + CS_vol + MB_vol) * f_meteo
        fr_THg = 4.0e-12 * 200.59 * 1e3
        fr_ratio = 0.04834
        fr_Hg2 = fr_THg * (1 - fr_ratio)
        fr_Hg0 = fr_THg * fr_ratio
        fr = np.array([fr_vol * fr_Hg0, fr_vol * fr_Hg2], dtype=float)

        H = np.exp(-2404.3 / (Tw + 273.15) + 6.92)
        n_fw = 4.2844e-5 + (0.157 * (Tw + 64.993) ** 2 - 91.296) ** (-1)
        n_sw_A = 1.541 + 1.998e-2 * Tw - 9.52e-5 * Tw**2
        n_sw_B = 7.974 - 7.561e-2 * Tw + 4.724e-4 * Tw**2
        n_sw = 10 * n_fw * (1 + n_sw_A * s_p / 1000 + n_sw_B * (s_p / 1000) ** 2)
        pi_w = 2.26
        M_w = 18.01
        V_b = 12.74
        D_sw_Hg = 7.4e-10 * (273.15 + Tw) * np.sqrt(pi_w * M_w) / n_sw / V_b**0.6
        v_sw = n_sw / p_sw
        Sc_Hg = v_sw / D_sw_Hg
        k_w = 0.25 * u_10**2 * (Sc_Hg / 660) ** (-0.5) * 0.01 / 3600

        Kp = 10**6.83
        C_spm = 21.4e-9
        F_Hg2 = 1 - (1 / (1 + Kp * C_spm))
        V_set = 24.9e-3 / 24 / 3600 * A

        UV_B_in = 0.00015787 * lat**2 - 0.026359 * lat + 1.1026
        UV_A_in = 0.01044 * lat**2 - 1.8153 * lat + 79.325
        PAR_in = (0.0217 * lat**2 - 3.8897 * lat + 178.7374) * 1e6 / 24 / 3600 / 4.6
        I_trans = 0.028
        UV_B = UV_B_in * (1 - ice) + UV_B_in * ice * I_trans
        UV_A = UV_A_in * (1 - ice) + UV_A_in * ice * I_trans
        PAR = PAR_in * (1 - ice) + PAR_in * ice * I_trans

        Chl = 0.0973 * s_p - 1.9444

        a_ph_443 = 10 ** (0.5877 * np.log10(Chl) - 1.5160)
        DOC = -3.85 * s_p + 190
        a_g_443 = (DOC - 55) / 357
        a_d_443 = 0.016
        a_tw_443 = 0.0145
        a_443 = a_ph_443 + a_g_443 + a_d_443 + a_tw_443

        Kd_490 = 0.0166 + 0.0733 * Chl**0.6715
        Kd_380 = 3.313 * a_443**2 + 0.887 * a_443 + 0.03
        Kd_340 = -1.965 * Kd_380**2 + 1.979 * Kd_380 - 0.009
        Kd_305 = -7.07 * Kd_380**2 + 4.424 * Kd_380 - 0.008
        Kd_PAR = 0.00665 + 0.874 * Kd_490 - 0.00121 * Kd_490**(-1)

        z_eu = np.log(100) / Kd_PAR
        D_irr = 16
        PAR_unit = PAR * 0.395
        Pb_opt = 2.0 * 10 ** (-8.7e-3 * (MLD / 2))

        NPP = 0.66125 * z_eu * Chl * D_irr * (PAR_unit / (PAR_unit + 4.1)) * Pb_opt
        pe_ratio = -0.0081 * Tw + 0.0668 * np.log(Chl / z_eu) + 0.426
        OCRR_unit = NPP / z_eu * (1 - pe_ratio)
        OCRR = OCRR_unit / 1000 / 24 / 12

        RHg_UVA = (-0.0274 * s_p**2 + 1.6344 * s_p - 23.752) / (1 + Kp * C_spm)
        RHg_UVB = 0.869 / (1 + Kp * C_spm)

        UV_B_exp = 1.4
        UV_A_exp = 6.2

        kr_UVA = (0.0088 * s_p**2 - 0.513 * s_p + 7.5412) / UV_A_exp
        kr_UVA_d = kr_UVA * UV_A * (1 - np.exp(-Kd_340 * MLD)) / Kd_340 / MLD
        k_UVA_2 = kr_UVA_d * RHg_UVA

        kr_UVB = 0.207 / UV_B_exp
        kr_UVB_d = kr_UVB * UV_B * (1 - np.exp(-Kd_305 * MLD)) / Kd_305 / MLD
        k_UVB_2 = kr_UVB_d * RHg_UVB

        if MLD >= z_eu:
            k_Bio_2 = 86 * OCRR
        else:
            k_Bio_2 = 86 * OCRR * MLD / z_eu

        k_2 = k_UVB_2 + k_UVA_2 + k_Bio_2 * 2

        UV_B_exp2 = 1.35
        UV_A_exp2 = 5.33

        ko_UVA = 0.340 / UV_A_exp2
        ko_UVB = 1.77 / UV_B_exp2

        ko_UVB_d = ko_UVB * UV_B * (1 - np.exp(-Kd_305 * MLD)) / Kd_305 / MLD
        ko_UVA_d = ko_UVA * UV_A * (1 - np.exp(-Kd_340 * MLD)) / Kd_340 / MLD
        k_UVB_1 = ko_UVB_d
        k_UVA_1 = ko_UVA_d

        k_Dark_1 = 1e-7 * 3600

        if MLD >= z_eu:
            k_Bio_1 = 140 * OCRR
        else:
            k_Bio_1 = 140 * OCRR * MLD / z_eu

        k_1 = k_UVB_1 + k_UVA_1 + k_Dark_1 + k_Bio_1
        k = np.array([k_1 / 3600, k_2 / 3600], dtype=float)

        diff_d = 50 / 2
        g = 9.8
        r_coeff = 0.2
        p_air = 1.3

        Cd_oa = 1.5e-3
        Cd_oi = 5e-3
        tou_oa = (1 - ice) * p_air * Cd_oa * u_10**2
        tou_oi = ice * p_sw * 1000 * Cd_oi * (0.015 * u_10) ** 2
        kappa = 0.41

        u_fric = np.sqrt((tou_oa + tou_oi) / (p_sw * 1000))
        e_diss = u_fric**3 / (kappa * MLD)

        N_bv = g / (p_sw * 1000) * (p_wml - p_sw * 1000) / diff_d
        Dw = r_coeff * e_diss * N_bv ** (-1)

        Cw = np.array([298e-15 * 200.59 * 1000, 826e-15 * 200.59 * 1000], dtype=float)

        c_ss = solve_cao_ss_direct(
            A,
            V,
            M_vol,
            MB,
            CS,
            ESS,
            CAO_nt,
            CAO_et,
            V_set,
            F_Hg2,
            Kp,
            k,
            H,
            k_w,
            ice,
            Cw,
            Dw,
            diff_d,
            V_A,
            L_ess,
            L_cs_v,
            L_w_v,
            L_n,
            L_e,
            L_v_v,
            V_dep,
            air_ratio,
            P,
            fr,
        )

        c1_ss = c_ss[0] * 1e9 / 200.59
        c2_ss = c_ss[1] * 1e9 / 200.59
        c3_ss = c_ss[2] * 1e9

        S0_HgIISink = V_set * c_ss[1] * (1 - F_Hg2) * Kp * 1e-6
        S0_Hg0toHgII = k[0] * V * c_ss[0]
        S0_HgIItoHg0 = k[1] * V * c_ss[1]
        AS = (1 - ice) * k_w * (c_ss[0] - c_ss[2] / H) * A
        Diff_Hg0 = A * Dw * (Cw[0] - c_ss[0]) / diff_d
        Diff_HgII = A * Dw * (Cw[1] - c_ss[1] * (1 - F_Hg2)) / diff_d
        Current_Hg0 = MB[0] + CS[0] + ESS[0]
        Current_HgII = MB[1] + CS[1] + ESS[1]

        NSE_num = (DGM_obs - c1_ss) ** 2
        NSE_den = (DGM_obs - 0.12056) ** 2
        PBIAS_num = (DGM_obs - c1_ss) * 100
        PBIAS_den = DGM_obs
        Net = S0_HgIItoHg0 - S0_Hg0toHgII

        rows.append(
            {
                "RowID": i + 1,
                "Temperature": Tw,
                "Salinity": s_p,
                "SeaIce": ice,
                "WindSpeed": u_10,
                "Hg0_water": c1_ss,
                "DGM_obs": DGM_obs,
                "HgII_water": c2_ss,
                "Hg0_air": c3_ss,
                "GEM_obs": GEM_obs,
                "MLD": MLD,
                "Dw": Dw,
                "P_sink": S0_HgIISink,
                "AS": AS,
                "Net": Net,
                "Diff_Hg0": Diff_Hg0,
                "Current_Hg0": Current_Hg0,
                "fr_Hg0": fr[0],
                "Diff_HgII": Diff_HgII,
                "Current_HgII": Current_HgII,
                "fr_HgII": fr[1],
                "Red_UVB": k_UVB_2,
                "Red_UVA": k_UVA_2,
                "Red_Bio": k_Bio_2,
                "Ox_UVB": k_UVB_1,
                "Ox_UVA": k_UVA_1,
                "Ox_Dark": k_Dark_1,
                "Ox_Bio": k_Bio_1,
                "a_443": a_443,
                "NSE_num": NSE_num,
                "NSE_den": NSE_den,
                "PBIAS_num": PBIAS_num,
                "PBIAS_den": PBIAS_den,
            }
        )

    result_table = pd.DataFrame(rows)
    result_table.to_csv(output_summary_csv, index=False)

    obs = result_table["DGM_obs"].to_numpy(dtype=float)
    sim = result_table["Hg0_water"].to_numpy(dtype=float)

    valid = np.isfinite(obs) & np.isfinite(sim)
    obs = obs[valid]
    sim = sim[valid]

    n = obs.size

    slope, intercept = np.polyfit(obs, sim, 1)
    sim_fit = slope * obs + intercept

    ss_res_fit = np.sum((sim - sim_fit) ** 2)
    ss_tot_fit = np.sum((sim - np.mean(sim)) ** 2)
    R2_reg = 1 - ss_res_fit / ss_tot_fit

    obs_mean = np.mean(obs)
    sim_mean = np.mean(sim)
    num = np.sum((obs - obs_mean) * (sim - sim_mean))
    den = np.sqrt(np.sum((obs - obs_mean) ** 2) * np.sum((sim - sim_mean) ** 2))

    r = num / den
    R2_corr = r**2

    RMSE = np.sqrt(np.mean((sim - obs) ** 2))
    PBIAS = 100 * np.sum(obs - sim) / np.sum(obs)
    Bias = np.mean(sim - obs)

    print("\n===== Model Performance Metrics =====")
    print(f"N             = {n:d}")
    print(f"R^2 (reg)     = {R2_reg:.4f}")
    print(f"R^2 (corr^2)  = {R2_corr:.4f}")
    print(f"r             = {r:.4f}")
    print(f"RMSE          = {RMSE:.4f}")
    print(f"Slope         = {slope:.4f}")
    print(f"Y-intercept   = {intercept:.4f}")
    print(f"PBIAS (%)     = {PBIAS:.2f}")
    print(f"Mean Bias     = {Bias:.4f}")

    metrics_table = pd.DataFrame(
        [
            {
                "N": n,
                "R2_reg": R2_reg,
                "R2_corr": R2_corr,
                "r": r,
                "RMSE": RMSE,
                "Slope": slope,
                "Intercept": intercept,
                "PBIAS_percent": PBIAS,
                "MeanBias": Bias,
            }
        ]
    )
    metrics_table.to_csv(output_metrics_csv, index=False)

    return result_table, metrics_table


if __name__ == "__main__":
    run_model()
