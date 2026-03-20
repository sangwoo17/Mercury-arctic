%HW1_ode.m

function c_dot = HW1_ode(t,c,A,D,k)


% Soil Hg(II) to Hg(0)
R_sred = k(1)*A(1)*(1-A(2))*A(3)*D(4)*c(1,3); %/s * m2  * m * kg/m3 = kg/s

% Litterfall of vegetation to soil
R_LII = k(2)*A(1)*(1-A(2))*A(3)*D(3)*c(1,4);   %/s * m2 * kg/m2 * kg/kg = kg/s 
R_L0 = k(2)*A(1)*(1-A(2))*A(3)*D(3)*c(2,4);    %/s * m2 * kg/m2 * kg/kg = kg/s   

%Redox in water
R_red = k(3)*A(1)*A(2)*D(2)*c(1,1);    %/s * m2 * m * kg/m3 = kg/s 
R_ox = k(4)*A(1)*A(2)*D(2)*c(2,1);     %/s * m2 * m * kg/m3 = kg/s 

%Mass transfer 
M1 = M_Hg0aw*A(1)*A(2)*c(2,1); %m/s * m2 * kg/m3 = kg/s -> Hg(0) air to water
M2 = M_HgIIaw*A(1)*A(2)*c(1,1); %m/s * m2 * kg/m3 = kg/s -> Hg(II) air to water
M3 = M_Hg0as*A(1)*(1-A(2))*A(3)*c(2,1); %m/s * m2 * kg/m3 = kg/s -> Hg(0) air to soil
M4 = M_HgIIas*A(1)*(1-A(2))*A(3)*c(1,1); %m/s * m2 * kg/m3 = kg/s -> Hg(II) air to soil 
M5 = M_Hg0HgIIav*A(1)*(1-A(2))*A(3)*c(2,1); %m/s * m2 * kg/m3 = kg/s -> Hg(0) air to vegetation
M6 = M_Hg0HgIIav*A(1)*(1-A(2))*A(3)*c(1,1); %m/s * m2 * kg/m3 = kg/s -> Hg(II) air to vegetation
M7 = M_Hg0wa*A(1)*A(2)*D(2)*c(2,2); %m/s * m2 * kg/m3 = kg/s -> Hg(0) water to air
M8 = M_Hg0sa*A(1)*(1-A(2))*A(3)*c(2,3); %m/s * m2 * kg/m3 = kg/s -> Hg(0) soil to air

%Runoff
R0_stow = V_r0*sqrt(A(1)*(1-A(2))*A(3))*D(4)*c(2,3); %m/s * m2 * kg/m3 = kg/s
RII_stow = V_rII*sqrt(A(1)*(1-A(2))*A(3))*D(4)*c(1,3); %m/s * m2 * kg/m3 = kg/s

%Uptake
R_up0 = V_up * A(1)*(1-A(2))*A(3)*c(1,3); %m/s * m2 * kg/m3 = kg/s
R_upII = V_up * A(1)*(1-A(2))*A(3)*c(2,3);	%m/s * m2 * kg/m3 = kg/s

%Air flow
R_in_air0 = U*C_in_air0*sqrt(A(1))*D(1)*10^12; %m/s * ng/m3 *m2 * kg/ng = kg/s
R_in_airII = U*C_in_airII*sqrt(A(1))*D(1)*10^12; %m/s * ng/m3 *m2 * kg/ng = kg/s

% Water flow 
R_in_water0 = Qwater*C_in_w0*10^-12; %m3/s * pg/m3 * kg/pg = kg/s
R_in_waterII = Qwater*C_inwII*10^-12; %m3/s * pg/m3 * kg/pg = kg/s
  


c_dot(1,1) = -M2 - M4 - M6 + R_in_airII;   % Hg(II) in air
c_dot(2,1) = -M1 -M3 - M5 + M7 + M8 +R_in_air0;   % Hg(0) in air 
c_dot(1,2) = R_ox - R_red + M2 + RII_stow + R_in_waterII;   % Hg(II) in water
c_dot(2,2) = R_red - R_ox + M1 - M7 + R0_stow + R_in_water0;   % Hg(0) in water
c_dot(1,3) = R_sred + R_LII + M4 - RII_stow - R_upII;   % Hg(II) in soil
c_dot(2,3) = -R_sred + R_L0 + M3 + M8 - R0_stow - R_up0;   % Hg(0) in soil
c_dot(1,4) = -R_LII + M6 + R_upII;   % Hg(II) in vegetation
c_dot(2,4) = -R_L0 + M5 + R_up0; % Hg(0) in vegetation

