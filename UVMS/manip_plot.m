function manip_plot(eta, q_sym, q_val, com, T, l, r, lt)

R_I_b = [ cos(eta(6))*cos(eta(5)) sin(eta(6))*cos(eta(5)) -sin(eta(5)); ...
    -sin(eta(6))*cos(eta(4)) + cos(eta(6))*sin(eta(5))*sin(eta(4)) cos(eta(6))*cos(eta(4)) + sin(eta(6))*sin(eta(5))*sin(eta(4)) sin(eta(4))*cos(eta(5));
    sin(eta(6))*sin(eta(4)) + cos(eta(6))*sin(eta(5))*cos(eta(4)) -cos(eta(6))*sin(eta(4)) + sin(eta(6))*sin(eta(5))*cos(eta(4)) cos(eta(4))*cos(eta(5)) ];

for k = 1:numel(q_sym)
    rot = double(subs(T{k}(1:3,1:3), q_sym, q_val));
    temp = eta;
    temp(1:3) = double(subs(com{k}, q_sym, q_val)) + temp(1:3);
    temp(4:end) = -[atan2(rot(3,2),rot(2,2)); -asin(rot(1,2)); atan2(rot(1,3),rot(1,1))] + eta(4:end);
    hold on
    cylindrical_plot(temp',l(k),r(k)*2,lt)
end 
hold off

end