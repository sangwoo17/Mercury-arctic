%Example2_ode.m

function c_dot = Example2_ode(t,c,V,H,A,Q_in,C_in,E,v_d,k_ox)
Q_out = Q_in;
F_out = Q_out * c;
D = A * (v_d/1000)*c; %mg/d, mass flux for deposition
L = k_ox * c * V; %mg/d, mass loss due to oxidation 
c_dot = (E + Q_in*C_in - D - L - F_out)/V;
