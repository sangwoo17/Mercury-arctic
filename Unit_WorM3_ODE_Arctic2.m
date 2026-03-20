%Unit_WorM3_ODE_v001.m

function c_dot = Unit_WorM3_ODE_Arctic2(t,c,A,V,D,M_vol,R,Co,L,CS,CAO_t,V_set,F_Hg2,k,H,k_w,Hg_air,ice,Cw,Dw,Diff_d)


%surface ocean 

Diff_Hg0 = A*Dw*(Cw(1)-c(1))/Diff_d;
Diff_HgII = A*Dw*(Cw(2)-c(2))/Diff_d;

M(1) = M_vol*c(1);
M(2) = M_vol*c(2);


AS = (1-ice)*k_w*(c(1)-Hg_air/H)*A;

CAO_Hg0 = CAO_t*c(1);
CAO_HgII = CAO_t*c(2);

S0_HgIISink = V_set*A*F_Hg2*c(2);

S0_Hg0toHgII = k(1)*V*c(1);
S0_HgIItoHg0 = k(2)*V*c(2);

c_dot(1,1) = (1/V)*(- S0_Hg0toHgII + S0_HgIItoHg0 + M(1) +R(1) +L(1)-CS(1)- CAO_Hg0 -AS + Diff_Hg0);
c_dot(2,1) = (1/V)*(- S0_HgIItoHg0 + S0_Hg0toHgII - S0_HgIISink + D(2)*(1-ice) + M(2) + R(2) + Co(2) + L(2)-CS(2) - CAO_HgII + Diff_HgII);