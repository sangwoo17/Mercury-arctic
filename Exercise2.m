%Exercise2

clc;

T = 0:5:30;
u10 = 0:0.1:20;

kw1 = 0.0025 * (u10.^2); 
kw2 = 0.0144 + 0.00144*(u10.^2);
plot(u10, kw1, 'LineWidth', 2);
hold on;
plot(u10, kw2, 'Linewidth',2);
hold on;

for i = 1:length(T)
    kw3 = 0.0025 * (u10.^2) * (0.0315 * T(i) + 0.784);
    plot(u10, kw3, 'LineWidth', 2);
         %'DisplayName', ['kw3, T = ' num2str(T(i)) '¡ÆC']);
end

hold off;
xlabel('u10 (m/s)');
ylabel('kw');
title('mass transfer velocity');
grid on;
xlim([0 30]);
legend('nitingale','schwarzenbach','0','5','10','15','20','25','30');
