close all
clear

%% ПАРАСЕТРЫ
L = 2500; % Длина ROV [mm]
H = 1500; % Высота ROV [mm]
W = 1600; % Ширина ROV [mm]

RotOX = @(angle)[1 0 0; 0 cos(angle) -sin(angle); 0 sin(angle) cos(angle)];
RotOY = @(angle)[cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle)];
RotOZ = @(angle)[cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];
% Rot = @(eta) ...
%         [cos(eta(6))*cos(eta(5)) sin(eta(6))*cos(eta(5)) -sin(eta(5)); ...
%         -sin(eta(6))*cos(eta(4)) + cos(eta(6))*sin(eta(5))*sin(eta(4)) ...
%         cos(eta(6))*cos(eta(4)) + sin(eta(6))*sin(eta(5))*sin(eta(4)) ...
%         sin(eta(4))*cos(eta(5));
%         sin(eta(6))*sin(eta(4)) + cos(eta(6))*sin(eta(5))*cos(eta(4)) ...
%         -cos(eta(6))*sin(eta(4)) + sin(eta(6))*sin(eta(5))*cos(eta(4)) ...
%         cos(eta(4))*cos(eta(5))];
Rot = @(eta)([RotOZ(eta(6))*RotOY(eta(5))*RotOX(eta(4))]);
eta = sym('eta', [6 1], 'real');

Tr_AUV = @(eta)[ [Rot(eta); zeros(1,3)] [eta(1:3); 1] ];




%% ПАРАМЕТРЫ МАНИПУЛЯТОРА
origin = [750; 0; -750].*10^-3; % координаты основания
n = 2; 
l = [1; 1]; 
r = [0.1; 0.1];
Rot0 = rotx(90); 

q = sym('q', [n 1], 'real'); 

c{1} = [-l(1)/2; 0; 0];
c{2} = [-l(2)/2; 0; 0];

%      a  alpha d   q
DH = [l(1)  0   0  q(1) - pi/2; ...
      l(2)  0   0  q(2)];     % DH Parameter Matrix

A = @(a, alpha, d, theta) ...
        ([cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha) a*cos(theta); ...
        sin(theta) cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta); ...
        0 sin(alpha) cos(alpha) d; ...
        0 0 0 1]);

% Homogeneous transformations solution
Tr = cell(n,1);
Tr0 = eye(4,'sym'); 
Tr0(1:3,1:3) = Rot0; 
Tr0(1:3,4) = origin;
Tr0 = Tr_AUV(eta) * Tr0;
T = Tr0;
for i = 1:n
    T = T * A(DH(i,1),DH(i,2),DH(i,3),DH(i,4));
    Tr{i} = simplify(T);
end

%% Центры масс
r_c_m = cell(n,1);
for i = 1:n    
    temp = Tr{i}*[[1;0;0;0] [0;1;0;0] [0;0;1;0] [c{i};1]];
	r_c_m{i} = temp(1:3,4);  
end

%% МОДЕЛИРОВАНИЕ ДВИЖЕНИЯ
eta_start = [0, 0, -5, 0, pi/4, pi/4]'; % начальная глубина 5 метров
eta_finish = [5, 5, -10, 0, pi/5, pi/7]'; 
q_start = deg2rad([30, 90]'); % начальные углы манипулятора
q_finish = deg2rad([-30, 90]');

%% ПОСТРОЕНИЯ ГРАФИКОВ
figure, title('Позиционирование AUV')
hold on

xlabel('x(t)'), ylabel('y(t)'), zlabel('z(t)'), grid on
set(gca, 'YDir', 'reverse');
view([-1,1,1])
axis equal
legend('off')

rectangular_plot(eta_start', L*10^-3, W*10^-3, H*10^-3, '--r')
xyz = double(subs(Tr0, [eta; q], [eta_start; q_start]));
x = xyz(1,4); y = xyz(2,4); z = xyz(3,4);
for i=1:numel(Tr)
    xyz = double(subs(Tr{i}, [eta; q], [eta_start; q_start]));
    x = [x; xyz(1,4)]; y = [y; xyz(2,4)]; z = [z; xyz(3,4)];
end
hold on
plot3(x,y,z,'r','LineWidth', 2.5)

rectangular_plot(eta_finish', L*10^-3, W*10^-3, H*10^-3, 'r')
xyz = double(subs(Tr0, [eta; q], [eta_finish; q_finish]));
x = xyz(1,4); y = xyz(2,4); z = xyz(3,4);
for i=1:numel(Tr)
    xyz = double(subs(Tr{i}, [eta; q], [eta_finish; q_finish]));
    x = [x; xyz(1,4)]; y = [y; xyz(2,4)]; z = [z; xyz(3,4)];
end
hold on
plot3(x,y,z,'r','LineWidth', 2.5)

% hold on
% manip_plot(eta(1,:)', q, q_vect(1,:)', com, Ti, l, r, '--r')
% hold on
% manip_plot(eta(end,:)', q, q_vect(end,:)', com, Ti, l, r, 'r')
