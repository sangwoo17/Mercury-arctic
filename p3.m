clear; clc;

k(1,1) = 0.5; % k1, h-1
k(2,1) = 2; % k-1, h-1
k(3,1) = 4; % k2, h-1
k(4,1) = 3; %k1', h-1
k(5,1) = 0;
k(6,1) = 0.5;


% Integration
tspan = [0 20]; % integration time span = 0 to 5000 days 
c0(1,1) = 280; % Hg(II)-ori
c0(2,1) = 0; % Hg(I)
c0(3,1) = 0; % Hg(0)
c0(4,1) = 0; % Hg(II)-mod

options = odeset('RelTol', 1e-13);
[t,c] = ode45(@p3_ode,tspan,c0,options,k);

figure
plot(t,c(:,3));
ylabel('Hg(0)');
xlabel ('Time (h)');


