function rectangular_plot(eta,L,W,H,lt)

dx = L/2;
dy = W/2;
dz = H/2;

R_I_b = [ cos(eta(6))*cos(eta(5)) sin(eta(6))*cos(eta(5)) -sin(eta(5)); ...
    -sin(eta(6))*cos(eta(4)) + cos(eta(6))*sin(eta(5))*sin(eta(4)) cos(eta(6))*cos(eta(4)) + sin(eta(6))*sin(eta(5))*sin(eta(4)) sin(eta(4))*cos(eta(5));
    sin(eta(6))*sin(eta(4)) + cos(eta(6))*sin(eta(5))*cos(eta(4)) -cos(eta(6))*sin(eta(4)) + sin(eta(6))*sin(eta(5))*cos(eta(4)) cos(eta(4))*cos(eta(5)) ];

A = eta(1:3) + [-dx dy -dz]; A1 = eta(1:3) + [-dx dy dz];
A = (A - eta(1:3)) * R_I_b + eta(1:3);
A1 = (A1 - eta(1:3)) * R_I_b + eta(1:3);
B = eta(1:3) + [-dx -dy -dz]; B1 = eta(1:3) + [-dx -dy dz];
B = (B - eta(1:3)) * R_I_b + eta(1:3); 
B1 = (B1 - eta(1:3)) * R_I_b + eta(1:3);
C = eta(1:3) + [dx -dy -dz]; C1 = eta(1:3) + [dx -dy dz];
C = (C - eta(1:3)) * R_I_b + eta(1:3); 
C1 = (C1 - eta(1:3)) * R_I_b + eta(1:3);
D = eta(1:3) + [dx dy -dz]; D1 = eta(1:3) + [dx dy dz];
D = (D - eta(1:3)) * R_I_b + eta(1:3); 
D1 = (D1 - eta(1:3)) * R_I_b + eta(1:3);

%plot3(eta(1), eta(2), eta(3), 'or')
%xlabel('X'), ylabel('Y'), zlabel('Z')
%set(gca, 'YDir', 'reverse');
%view([-1,1,1]), axis equal, grid on, 
hold on
plot3([A(1) B(1) C(1) D(1) A(1)], ... 
    [A(2) B(2) C(2) D(2) A(2)], ...
    [A(3) B(3) C(3) D(3) A(3)], lt)
plot3([A1(1) B1(1) C1(1) D1(1) A1(1)], ...
    [A1(2) B1(2) C1(2) D1(2) A1(2)], ...
    [A1(3) B1(3) C1(3) D1(3) A1(3)], lt)
plot3([A(1) A1(1)], [A(2) A1(2)], [A(3) A1(3)], lt)
plot3([B(1) B1(1)], [B(2) B1(2)], [B(3) B1(3)], lt)
plot3([C(1) C1(1)], [C(2) C1(2)], [C(3) C1(3)], lt)
plot3([D(1) D1(1)], [D(2) D1(2)], [D(3) D1(3)], lt)
hold off

end

