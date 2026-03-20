%Example3_ode.m

function c_dot = Example3_ode(t,c,V,Q,C_in,E)

% air
c_dot(1,1) = ( Q(1,1)*C_in(1,1) + E(1,1) - Q(1,1)*c(1,1) ) / V(1,1) ;


% water
c_dot(2,1) = ( Q(2,1)*C_in(2,1) + E(2,1) - Q(2,1)*c(2,1) ) / V(2,1) ;