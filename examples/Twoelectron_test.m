clear
clc

%% 실험 데이터 (예시)
t_exp = [0 0.5 1 2 3 4 5]; % 시간
HgI_exp = [0 20 35 50 45 30 20]; % Hg(I) 농도 데이터 예시

%% 초기 농도
HgII0 = 100; % 예: 100 nM
HgI0 = 0;
Hg0_0 = 0;
c0 = [HgII0; HgI0; Hg0_0];

%% ODE 함수 핸들
odefun = @(t, c, k) twoelectron_ode(t, c, k);

%% Objective Function (Sum of Squared Errors)
obj_fun = @(k) sum_of_squares(k, t_exp, HgI_exp, c0, odefun);

%% 초기값 추정
k0 = [1, 0.5, 0.1]; % k1, k2, k3 초기값

%% Optimization (fminsearch)
k_fit = fminsearch(obj_fun, k0);

%% 결과 확인
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



