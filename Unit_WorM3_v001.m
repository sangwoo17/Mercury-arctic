clc; clear; 

%Unit_WorM_v001.m

A(1) = 510e12; %m2 (1) = lower atmosphere
A(2) = A(1)*0.7; %m2 
A(3) = A(1)-A(2);

D(1) = 6000;
D(2) = 100;
D(3) = 0.15;

V(1) = A(1)*D(1);
V(2) = A(2)*D(2);
V(3) = (A(1)-A(2))*D(3);
CR02 = 100; %Concentration ratio of gaseous Hg(0) and gaseous Hg(II) in the atmosphere

Tw = 282; %temperature of water in K
H = exp(-4633/Tw + 14.53); % Dimensinless H (conc-gas/con-aq) (Waenberg et al., 2001)
%Kaw = H;

%DEPOSITION VELOCITIES

Vd(1,1) = 2.1e-5; %m/s
Vd(2,1) = 0.011; %m/s
Vd(3,1) = 4.62e-5; %m/s Hg(0) air to soil
Vd(3,2) = 0.005; %m/s Hg(II) air to soil

Vd(3,2) = Vd(3,2) * 1.001;

Vsink = 6.1e-12; %2.2e-8 mh-1; %6.1e-12 ms-1 % from The tool - ocean particle sinking

%RATE CONSNTANTS FOR REACTIONS IN WATER ONLY

k(1) = 0.1; % (/h) rate constant for Hg(II) reduction
k(2) = 1.9; % (/h) rate constant for Hg(0) oxidation
k(3) = 1e-10; %/s - rate constant for HgII reduction in soil 


k(1) = k(1)/D(2)*(1-exp(-100));
k(2) = k(2)/D(2)*(1-exp(-100));


k(1) = k(1)/(60*60); %/s
k(2) = k(2)/(60*60); %/s




f_p = 2e-6; % volume fraction of suspended solids in water
Ksus_w(1) = 3e4; % partition coefficient for Hg(0) in suspended solids in water
Ksus_w(2) = 8e4; % partition coefficient for Hg(II) in suspended solids in water

%pre-industrial ~ 9000 years

E = 64.3/3; %Emissions of Hg(0) + Hg(II) in gs-1 to lower atmosphere = 2000 t/y, Pacyna 2000 in

% c(1) = Hg0 in air; c(2) = HgII in air; c(3) = bulk Hg0 in water, c(4) = %
% bulk HgII in water c(5) = Hg(0) in soil c(6) = Hg(II) in soil
C0 = [2e-15 1e-15 1e-15 1e-15 1e-15 1e-15];

options = odeset('RelTol', 1e-13, 'AbsTol', [2e-18 1e-18 1e-18 1e-18 1e-18 1e-18], 'NonNegative', [1 2 3 4 5 6]);

tmax = 9000 * 365 * 24 * 60 * 60; % 9000 years
tspan = [0 tmax];

[PreIndT, PreIndConc] = ode15s(@Unit_WorM3_ODE_v001,tspan,C0,options,A,D,V,CR02,H,Vd,Vsink,k,f_p,Ksus_w,E);

figure
plot(PreIndT/(365*24*60*60), PreIndConc(:,1),'-o')
xlabel('Pre-industrial time (years)')
ylabel('Hg(0) concentration in air (g/m3)')


ElementCount_PreIndConc = length(PreIndConc);

%post-industrial ~ 50 years

E = 64.3; % Emissions of Hg(0) + Hg(II) in gs-1 to lower atmosphere = 2000t/y, Pacyna, 2000
C0(1,:) = PreIndConc(ElementCount_PreIndConc,:);


options = odeset('RelTol', 1e-13,'AbsTol', [1e-18 1e-18 1e-18 1e-18 1e-18 1e-18], 'NonNegative', [1 2 3 4 5 6]);

tmax = 50*365*24*60*60; % 50years
tspan = [0 tmax];

[T,C] = ode15s(@Unit_WorM3_ODE_v001,tspan,C0,options,A,D,V,CR02,H,Vd,Vsink,k,f_p,Ksus_w,E);

figure
plot(T/(365*24*60*60),C(:,1),'-o')
xlabel('Post-Indsutrial time (years)')
ylabel('Hg(0) concentration in air (g/m3)')
