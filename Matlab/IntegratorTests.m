function IntegratorTests

% Can delete later, test of RK3CN variable time-step integrator 
errstore = [];
tvec = [0.0001 0.001 0.002 0.005 0.01 ];%0.02 0.05 0.1 0.2];
for dt = tvec
y = [0;2;20];
% y = [0;2];
t_final = 2*pi;
t=0;
ystore = [t;y];

a0=8/15;a1=5/12;a2=3/4;
b1=-17/60;b2=-5/12;

f = @lorenz;
while t<t_final
    
%     y = y + f(y,t)*dt;
    rhs = f(y,t);
    y = y + dt*a0*rhs;
    rhs2 = f(y,t);
    y = y + dt*(a1*rhs2+ b1*rhs);
    rhs = f(y,t);
    y = y + dt*(a2*rhs + b2*rhs2);
    
    t=t+dt;
    
    ystore = [ystore [t;y]];
end

% plot3(ystore(1,:),ystore(2,:),ystore(3,:))
% plot(ystore(1,:),ystore(2,:))
if dt == 0.0001
    solstore = ystore;
else
    ytrue = solstore(2:end,1:dt/0.0001:end);
    errstore = [errstore sum(sum(sqrt((ytrue-ystore(2:end,1:end-1)).^2)))/length(ystore)];
end



end
% Checked and this is 3rd order in time
plot(log10(tvec(2:end)),log10(errstore))

end

function out = lorenz(y,t)

out = [10*(y(2)-y(1)); y(1)*(28-y(3))-y(2); y(1)*y(2)-8/3*y(3)];
end


function out = HO(y,t)

out = [y(2); -y(1)];
end
