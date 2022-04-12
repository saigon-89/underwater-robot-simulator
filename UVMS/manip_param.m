origin = [750; 0; -750].*10^-3;
n = 2; % число звеньев
l = [1; 1]; % длины звеньев
r = [0.1; 0.1]; % радиусы звеньев
m = [10; 10]; % mass of each link [kg]
R0 = rotx(90); % начальная матрица поворота

q = sym('q', [n 1], 'real'); % обобщенные координаты (углы соединений)
g_accel = 9.81;

% Координаты центров масс каждого звена в его собственной системе отсчета
c{1} = [-l(1)/2; 0; 0];
c{2} = [-l(2)/2; 0; 0];

%      a  alpha d   q
DH = [l(1)  0   0  q(1)-pi/2; ...
      l(2)  0   0  q(2)];     % DH Parameter Matrix

% https://automaticaddison.com/how-to-find-denavit-hartenberg-parameter-tables/
% DH = [0    pi/2  l(1)  q(1); ...
%      l(2)   0     0    q(2)];     % DH Parameter Matrix

% inertia tensor for each link relative to the inertial frame stored in an nx1 cell array
I = cell(1,n);
I{1} = m(1).*diag([0.5*r(1)^2, (3*r(1)^2 + l(1)^2)/12, (3*r(1)^2 + l(1)^2)/12]);
I{2} = m(2).*diag([0.5*r(2)^2, (3*r(2)^2 + l(2)^2)/12, (3*r(2)^2 + l(2)^2)/12]);