function ssq = sum_of_squares(k, t_exp, HgI_exp, c0, odefun)
    [t_sol, c_sol] = ode45(@(t, c) odefun(t, c, k), t_exp, c0);
    HgI_sim = c_sol(:, 2);
    ssq = sum((HgI_sim - HgI_exp(:)).^2);
end
