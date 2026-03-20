% Example2.m

% Environment

V = 51; % m3
H = 3; % m
A = V/H; % m2

% Sources
    % Emissions
    E = 1.5; % mg/d
    
    % Inflow
    Q_in = 100; % m3/d
    C_in = 0.01; % mg/m3
    F_in = Q_in * C_in; %mg/d

% Sinks 

    % Deposition
    
    v_d = 1; % mm/d
    v_d = v_d/1000; %m/d
  
    
    % Oxidation loss 
    t_ox = 20; %d
    k_ox = 1/t_ox; %1/d
   
    
    % Outflow 
    Q_out = Q_in;
    

    
% Integration
tspan = [0 5]; % integration time span = 0 to 5 days 
c0 = 1e-10; % initial concentration in the room (low C0 value is not equal to zero)

options = odeset('RelTol', 1e-3);
[t,c] = ode45(@Example2_ode,tspan,c0,options,V,H,A,Q_in,C_in,E,v_d,k_ox);

plot(t,c)

%Writing values
filename = 'Example2_results.xlsx';
writematrix(t,filename,'sheet',1,'Range','A2');
writematrix(c,filename,'sheet',1,'Range','B2');



    
