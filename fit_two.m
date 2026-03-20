clear
clc

%% === 실험 데이터 입력 ===

data = readmatrix('HgI_data.xlsx');
t_exp = data(:, 1);       % 시간 (h)
HgI_exp = data(:, 2);    % 실험 Hg(I) 농도 (nM)

%% === 초기 농도 조건 ===
HgII0 = 112;  % 초기 Hg(II) 농도
HgI0 = 0;
Hg0_0 = 0;
c0 = [HgII0; HgI0; Hg0_0];

%% === ODE 함수 핸들 ===
odefun = @(t, c, k) two_ode(t, c, k);

%% === Objective Function ===
obj_fun = @(k) sum_of_squares(k, t_exp, HgI_exp, c0, odefun);

%% === k1, k2, k3 초기값 ===
k0 = [1, 1, 0.5];

%% === Optimization (fminsearch) ===
k_fit = fminsearch(obj_fun, k0);

%% === ODE 해석 후 결과 확인 ===
[t_sol, c_sol] = ode45(@(t, c) odefun(t, c, k_fit), t_exp, c0);
HgI_sim = c_sol(:, 2);

figure;
plot(t_exp, HgI_exp, 'ro', 'DisplayName', 'Exp [Hg(I)]');
hold on
plot(t_exp, HgI_sim, 'b-', 'DisplayName', 'Fitted [Hg(I)]');
xlabel('Time (h)')
ylabel('[Hg(I)] Concentration (nM)')
legend show
title(sprintf('Fitted k1 = %.3f, k2 = %.3f, k3 = %.3f', k_fit(1), k_fit(2), k_fit(3)))
