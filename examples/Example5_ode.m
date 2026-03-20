%Example5_ode.m

function c_dot = Example5_ode(t,c,k)


c_dot(1,1) = - k(1,1)*c(1,1);
c_dot(2,1) = - k(2,1)*c(2,1) + k(1,1)*c(1,1) + k(3,1)*c(3,1);
c_dot(3,1) = - k(3,1)*c(3,1) + k(2,1)*c(2,1);