%Unit_WorM3_ODE_v001.m

function c_dot = Unit_WorM3_ODE_CAO(t,c,A,V,D,M_vol,MB,CS,ESS,CAO_nt,CAO_et,V_set,F_Hg2,k,H,k_w,Hg_air,ice,Cw,Dw,Diff_d)


%surface ocean 

M(1) = M_vol*c(1);
M(2) = M_vol*c(2);

Diff_Hg0 = A*Dw*(Cw(1)-c(1))/Diff_d;
Diff_HgII = A*Dw*(Cw(2)-c(2)*(1-F_Hg2))/Diff_d;

AS = (1-ice)*k_w*(c(1)-Hg_air/H)*A;

CAOn(1) = CAO_nt*c(1);
CAOn(2) = CAO_nt*c(2);

CAOe(1) = CAO_et*c(1);
CAOe(2) = CAO_et*c(2);

S0_HgIISink = V_set*A*F_Hg2*c(2);

S0_Hg0toHgII = k(1)*V*c(1);
S0_HgIItoHg0 = k(2)*V*c(2);


c_dot(1,1) = (1/V)*(- S0_Hg0toHgII + S0_HgIItoHg0 + M(1)  + MB(1)+CS(1)+ESS(1)- CAOn(1) - CAOe(1) -AS + Diff_Hg0);
c_dot(2,1) = (1/V)*(- S0_HgIItoHg0 + S0_Hg0toHgII - S0_HgIISink + D(2) + M(2)  + MB(2)+CS(2)+ESS(2) - CAOn(2) - CAOe(2) + Diff_HgII);