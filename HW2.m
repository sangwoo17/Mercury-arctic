A(1) = 510e12; %m2 (1) = lower atmosphere
A(2) = A(1)*0.7; %m2 

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

Vsink = 6.1e-12; %2.2e-8 mh-1; %6.1e-12 ms-1 % from The tool - ocean particle sinking

%RATE CONSNTANTS FOR REACTIONS IN WATER ONLY

k(1) = 0.15; % (/h) rate constant for Hg(II) reduction
k(2) = 1.9; % (/h) rate constant for Hg(0) oxidation

k(1) = k(1)/(60*60); %/s
k(2) = k(2)/(60*60); %/s

k(1) = k(1)/D(2)*(1-exp(-100));
k(2) = k(2)/D(2)*(1-exp(-100));

f_p = 2e-6; % volume fraction of suspended solids in water
Ksus_w(1) = 3e4; % partition coefficient for Hg(0) in suspended solids in water
Ksus_w(2) = 8e4; % partition coefficient for Hg(II) in suspended solids in water

c(1) = 1.0231e-9;
c(2) = 0;
c(3) = 9.8937e-09;
c(4) = 1.3865e-07;
c(5) = 1.0231e-09;
c(6) = 1.33e-10;

Hg0_AtoS0=Vd(1,1)*A(2)*(c(1)/H);

Hg0_S0toA=Vd(1,1)*A(2)*c(3)/((1-f_p)+f_p*Ksus_w(1)); %

HgII_AtoS0=Vd(2,1)*A(2)*c(1)/CR02;

S0_Hg0Sink = Vsink*A(2)*Ksus_w(1)*c(3)/((1-f_p)+f_p*Ksus_w(1));
S0_HgIISink = Vsink*A(2)*Ksus_w(2)*c(4)/((1-f_p)+f_p*Ksus_w(2));
S0_Hg0toHgII = k(2)*V(2)*(1-f_p)*c(3)/((1-f_p)+f_p*Ksus_w(1));
S0_HgIItoHg0 = k(1)*V(2)*(1-f_p)*c(4)/((1-f_p)+f_p*Ksus_w(2));


Hg0_AtoS = Vd(3,1)*(A(1)-A(2))*c(1);
HgII_AtoS = Vd(3,2)*(A(1)-A(2))*c(2);
Hg0_StoA=Vd(1,1)*A(3)*c(5);

R = (E-S0_Hg0Sink-S0_HgIISink)/E;