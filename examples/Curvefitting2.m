% ultra_minimal_fit_plot.m
clear; A0=201.93; c0=[A0;0;0]; s=1; % s=1/(5/60) if data=pg/5min
T=readtable('test.csv'); t=T.time; y=T.hg;
theta0=log([.5;1;.5]); theta=fminsearch(@(th)obj(th,t,y,c0),theta0);
k=exp(theta)

[~,c]=ode45(@(tt,cc)f(cc,k),t,c0);
yhat=k(2)*c(:,2);
R2=1-sum((y-yhat).^2)/sum((y-mean(y)).^2);
fprintf('R©÷=%.4f\n',R2);

plot(t,y,'o',t,yhat,'-'); xlabel('Time (h)'); ylabel('d[Hg(0)]/dt');
legend('data','fit'); title(sprintf('Fit result (R©÷=%.3f)',R2));

function S=obj(th,t,y,c0)
k=exp(th(:)); [~,c]=ode45(@(tt,cc)f(cc,k),t,c0);
S=sum((k(2)*c(:,2)-y).^2);
end

function c_dot=f(c,k)
c_dot=[-k(1)*c(1)+k(3)*c(2); k(1)*c(1)-(k(2)+k(3))*c(2); k(2)*c(2)];
end
