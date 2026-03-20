% Example5.m

clear;

k(1,1) = 0.5; % k1', h-1
k(2,1) = 1.5; % k2', h-1
k(3,1) = 0.1; % k3', h-1


% Integration
tspan = [0 20]; % integration time span = 0 to 5000 days 
c0(1,1) = 200; % Hg(II)
c0(2,1) = 0; % Hg(0)
c0(3,1) = 0; % Hg*
options = odeset('RelTol', 1e-3);
[t,c] = ode45(@Example5_ode,tspan,c0,options,k);

% concentrations in air
figure
plot(t,c(:,1));
ylabel('Hg(II)');
xlabel ('Time (h)');

figure
plot(t,c(:,2));
ylabel('Hg(0)');
xlabel ('Time (h)');


figure
plot(t,c(:,3));
ylabel('Hg*');
xlabel ('Time (h)');