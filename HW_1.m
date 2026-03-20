% HW_1.m

clear;

A(1) = 200e6;		% Total Environmental Area (m2) ? assume a square shape
A(2) = 0.7;		% Percent area under the lake
A(3) = 0.5;		% Percent vegetation area occupying soil

D(1) = 1000;    % Atmosphere Height (m)
D(2) = 20;		% Water depth (m)
D(3) = 0.2;		% Vegetation density (kg/m2)
D(4) = 0.2;		% Soil depth (m) 

k(1) = 1e-10;	% Rate constant (/s) Soil Hg(II) to Hg(0)
k(2) = 3.17e-8;	% Rate constant (/s) Litterfall of vegetation to soil, for Hg(0)/Hg(II)
k(3) = 1.0*1/3600;         % Rate constant (/s) Hg(0) to Hg(II)
k(4) = 0.5*1/3600;         % Rate constant (/s) Hg(II) to Hg(0)

RT = 10 * 365 * 24 *60 * 60;         % Residence time of water in lake, y * s/y = s
Qwater = A(1)*A(2)*D(2)/RT;           % Water flow rate, m3/s

% Mass transfer velocities (cm/s)

M_Hg0aw = 2e-3*0.01;       %Hg(0) air to water m/s
M_HgIIaw = 1.2*0.01;        %Hg(II) air to water m/s
M_Hg0as = 0.004*0.01;       %Hg(0) air to soil m/s
M_HgIIas = 0.6*0.01;        %Hg(II) air to soil m/s
M_Hg0HgIIav = 0.4*0.01;   %Hg(0)/Hg(II) air to vegetation m/s
M_Hg0wa = 0.4*0.01;         %Hg(0) water to air m/s
M_Hg0sa = 0.01*0.01;        %Hg(0) soil to air 	m/s

V_r0 = 1.08e-08;       % Runoff from soil to water, Hg(0) (m/s)
V_rII = 6.39e-12;      % Runoff from soil to water, Hg(II) (m/s)
V_up = 2.22e-07;       % Uptake of water by vegetation from soil, for both Hg(0)/ Hg(II) (m/s)
U = 5;                  	% Wind velocity (in and out) (m/s)
C_in_air0 = 2;              % Hg(0) in incoming air (ng/m3)
C_in_airII = 0.1;           % Hg(II) in incoming air (ng/m3)
C_in_w0 = 2;                % Hg(0) in incoming water (pg/m3)
C_in_wII = 5;           	% Hg(II) in incoming water (pg/m3)


% Integration
tspan = [0 20]; % integration time span = 0 to 20 seconds

c0(1,1) = 200;  % Hg(II) in air
c0(2,1) = 0;    % Hg(0) in air 
c0(1,2) = 0;    % Hg(II) in water
c0(2,2) = 0;    % Hg(0) in water
c0(1,3) = 0;    % Hg(II) in soil
c0(2,3) = 0;    % Hg(0) in soil
c0(1,4) = 0;    % Hg(II) in vegetation
c0(2,4) = 0;    % Hg(0) in vegetation

options = odeset('RelTol', 1e-3);
[t,c] = ode45(@HW1_ode,tspan,c0,options,A,D,k);



