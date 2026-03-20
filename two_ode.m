function c_dot = two_ode(t, c, k)
    c_dot = zeros(3, 1);
    c_dot(1,1) = -k(1)*c(1) + k(3)*c(2);    % Hg(II)
    c_dot(2,1) =  k(1)*c(1) - (k(2) + k(3))*c(2); % Hg(I)
    c_dot(3,1) =  k(2)*c(2);                % Hg(0)
end
