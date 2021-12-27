function dy = odefcn(t,y,M,C,D,g,tau)
    dy = zeros(12,1);
    dy(1:6) = y(7:end);
    dy(7:end) = - M\C(y(7:end))*y(7:end) - M\D(y(7:end))*y(7:end) - ...
        M\g(y(1:6)) + M\tau(t);
end