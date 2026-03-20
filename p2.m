clear;

k(1,1) = 1; % k1, h-1
k(2,1) = 12; % k-1, h-1
k(3,1) = 12; % k2, h-1
k(4,1) = 0.2; %k1', h-1
% Integration
tspan = [0 20]; % integration time span = 0 to 5000 days 
c0(1,1) = 195.24; % Hg(II)-ori
c0(4,1) = 0; % Hg(II)-mod
c0(2,1) = 0; % Hg(I)
c0(3,1) = 0; % Hg(0)
options = odeset('RelTol', 1e-13);
[t,c] = ode45(@p2_ode,tspan,c0,options,k);

%figure
%plot(t,c(:,2));
%ylabel('Hg(I)');
%xlabel ('Time (h)');
%hold on; 

figure
plot(t,c(:,2));
ylabel('Hg(0)');
xlabel ('Time (h)');


