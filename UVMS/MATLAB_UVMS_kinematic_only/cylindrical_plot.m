function cylindrical_plot(eta,L,D,lt)

dx = L/2;
r = D/2;

plot3(eta(1), eta(2), eta(3), 'or')
xlabel('X'), ylabel('Y'), zlabel('Z')
set(gca, 'YDir', 'reverse');
view([-1,1,1]), axis equal, grid on, 

R_I_b = [ cos(eta(6))*cos(eta(5)) sin(eta(6))*cos(eta(5)) -sin(eta(5)); ...
    -sin(eta(6))*cos(eta(4)) + cos(eta(6))*sin(eta(5))*sin(eta(4)) cos(eta(6))*cos(eta(4)) + sin(eta(6))*sin(eta(5))*sin(eta(4)) sin(eta(4))*cos(eta(5));
    sin(eta(6))*sin(eta(4)) + cos(eta(6))*sin(eta(5))*cos(eta(4)) -cos(eta(6))*sin(eta(4)) + sin(eta(6))*sin(eta(5))*cos(eta(4)) cos(eta(4))*cos(eta(5)) ];

theta = -pi:0.01:pi;

c = [zeros(1,numel(theta)) - dx; r*cos(theta); r*sin(theta)]';
c1 = [zeros(1,numel(theta)) + dx; r*cos(theta); r*sin(theta)]';

for i=1:length(theta)   
    c(i,:) = c(i,:) * R_I_b + eta(1:3); 
    c1(i,:) = c1(i,:) * R_I_b + eta(1:3); 
end

hold on
plot3(c(:,1),c(:,2),c(:,3), lt)
plot3(c1(:,1),c1(:,2),c1(:,3), lt)

for i=1:37:length(theta)
    plot3([c(i,1) c1(i,1)],[c(i,2) c1(i,2)],[c(i,3) c1(i,3)], lt)
end
hold off

end