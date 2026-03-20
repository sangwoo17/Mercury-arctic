% Exercise 1
clc;

c = 0:0.1:30;
T = 273.15 + c; % Temperature

K= exp(-2404.3 ./ T + 6.92); % Andersson 2008
H= exp(-4633 ./ T + 14.53); % Wangberg 2001

RPD = (K-H)./K*100;  % RPD

plot(c, K, 'LineWidth', 2);
hold on;
plot(c, H, 'Linewidth',2);
hold off;
xlabel('Temperature (¡ÆC)');
ylabel('K_{AW}');
title('Variation of K_{AW} with Temperature');
grid on;
xlim([0 30]);



plot(c, RPD);
xlabel('c');
ylabel('RPD');
title('RPD');
grid on;
xlim([0 30]);

