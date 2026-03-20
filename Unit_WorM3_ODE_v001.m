%Unit_WorM3_ODE_v001.m

function c_dot = Unit_WorM3_ODE_v001(t,c,A,D,V,CR02,H,Vd,Vsink,k,f_p,Ksus_w, E)

%Air

Hg0_AtoS0=Vd(1,1)*A(2)*(c(1)/H);

Hg0_S0toA=Vd(1,1)*A(2)*c(3)/((1-f_p)+f_p*Ksus_w(1)); %

HgII_AtoS0=Vd(2,1)*A(2)*c(1)/CR02;



%surface ocean 

S0_Hg0Sink = Vsink*A(2)*Ksus_w(1)*c(3)/((1-f_p)+f_p*Ksus_w(1));
S0_HgIISink = Vsink*A(2)*Ksus_w(2)*c(4)/((1-f_p)+f_p*Ksus_w(2));
S0_Hg0toHgII = k(2)*V(2)*(1-f_p)*c(3)/((1-f_p)+f_p*Ksus_w(1));
S0_HgIItoHg0 = k(1)*V(2)*(1-f_p)*c(4)/((1-f_p)+f_p*Ksus_w(2));

%soil

Hg0_AtoS = Vd(3,1)*(A(1)-A(2))*c(1);
HgII_AtoS = Vd(3,2)*(A(1)-A(2))*c(1)/CR02;
Hg0_StoA=Vd(3,1)*A(3)*c(5);
S_HgIItoHg0 = k(3)*c(6)*A(3)*D(3);

c_dot(1,1) = (1/V(1))*(-Hg0_AtoS0 + Hg0_S0toA - HgII_AtoS0 + E - Hg0_AtoS + Hg0_StoA - HgII_AtoS);
c_dot(2,1) = 0;
c_dot(3,1) = (1/V(2))*(Hg0_AtoS0 - Hg0_S0toA - S0_Hg0toHgII + S0_HgIItoHg0 - S0_Hg0Sink);
c_dot(4,1) = (1/V(2))*(HgII_AtoS0 - S0_HgIItoHg0 + S0_Hg0toHgII - S0_HgIISink);
c_dot(5,1) = (1/V(3))*(Hg0_AtoS-Hg0_StoA+S_HgIItoHg0);
c_dot(6,1) = (1/V(3))*(HgII_AtoS-S_HgIItoHg0);
