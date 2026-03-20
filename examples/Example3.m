% Example3.m

M = 200; %g/mol
%1: air 
%2: water
%3: sediment

% Environment
V(1,1) = 200000; % m3
V(2,1) = 1500; % m3
V(3,1) = 15; % m3

% Inflow
Q(1,1) = 1000; %m3/h (air inflow and outflow)
Q(2,1) = 1; %m3/h (water inflow and outflow)

C_in(1,1) = 0.01; % g/m3
C_in(2,1) = 0.1; % g/m3

% Sources
E(1,1) = 1; % mol/h - Emissions to air
E(1,1) = E(1,1) * M; % Emissions to air in g/h

E(2,1) = 0.24; % mol/h - Emissions to water
E(2,1) = E(2,1) * M; % Emissions to water in g/h

% Integration
tspan = [0 5000]; % integration time span = 0 to 5000 days 
c0(1,1) = 1e-10; % initial concentration in air
c0(2,1) = 1e-10; % initial concentration in water
options = odeset('RelTol', 1e-3);
[t,c] = ode45(@Example3_ode,tspan,c0,options,V,Q,C_in,E);

% concentrations in air
figure
plot(t,c(:,1));
ylabel('Concentration in air (g/m3)');
xlabel ('Time (h)');

% concentrations in water
figure
plot(t,c(:,2));
ylabel('Concentration in water (g/m3)');
xlabel ('Time (h)');

