% close all
clear

%% РАСЧЕТ ДИНАМИКИ ПОДВОДНОГО АППАРАТА %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ЗАГРУЗИТЬ ПАРАМЕТРЫ ИЗ ФАЙЛА
rov_param % [6]

%% ХАРАКТЕРНО ДЛЯ ПРЯМОУГОЛЬНОЙ ПРИЗМЫ
I0  = diag(m*[(W^2+H^2) (L^2+H^2) (W^2+L^2)]/12).*10^-6; % тензор инерции (1.6)

%% ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
% преобразование в кососимметричную матрицу
S = @(x)([ 0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0 ]); % (1.5)
r_g_c = r_g_c / 10^3; % переводим [mm] в [m]
r_b_c = r_b_c / 10^3; % переводим [mm] в [m]

%% РАСЧЕТ МАТРИЦЫ M
% расчет M_RB
M_RB = [ m*eye(3) -m*S(r_b_c); m*S(r_b_c) I0 ]; % (1.4)
% расчет M_A
M_A = rectangular_added_mass(L, H, W, rho, PF, PS, PT); % [6]
% расчет M
M = M_RB + M_A; % (1.3)

%% РАСЧЕТ C(v)
% расчет C_RB(v)
C_RB = @(v)([ zeros(3) -m*S(v(1:3)); -m*S(v(1:3)) -S(diag(I0).*v(4:end)) ]); % (1.7)
% расчет C_A(v)
diag_M_A = diag(M_A);
C_A = @(v)([ zeros(3) -S(diag_M_A(1:3).*v(1:3)); % (1.8)
    -S(diag_M_A(1:3).*v(1:3)) -S(diag_M_A(4:end).*v(4:end))]);
% расчет C(v)
C = @(v)(C_RB(v) + C_A(v)); % (1.3)

%% РАСЧЕТ g(n)
B = m * 9.81; % (1.9)
g = @(eta)([ 0; 
    0; 
    0;
    -(r_g_c(2)*B-r_b_c(2)*B)*cos(eta(5))*cos(eta(4)) + (r_g_c(3)*B-r_b_c(3)*B)*cos(eta(5))*sin(eta(4));
    (r_g_c(3)*B-r_b_c(3)*B)*sin(eta(5)) + (r_g_c(1)*B-r_b_c(1)*B)*cos(eta(5))*cos(eta(4));
    -(r_g_c(1)*B-r_b_c(1)*B)*cos(eta(5))*sin(eta(4)) - (r_g_c(2)*B-r_b_c(2)*B)*sin(eta(5)) ]); % (1.10)

%% РАСЧЕТ D(v)
[D_LIN, D_QUAD] = rectangular_damping(L, H, W, rho, PF, PS, PT, M_RB, M_A, B, r_g_c, r_b_c);
D = @(v)(D_LIN + D_QUAD.*abs(v)); % [6]

%% ВЫВОД РЕЗУЛЬТАТОВ
disp('Матрица M:')
disp(vpa(M, 4))

syms u v w p q r
disp('Матрица С(v):')
disp(vpa(C([u; v; w; p; q; r]), 4))

disp('Матрица D(v):')
disp(vpa(D([u; v; w; p; q; r]), 4))

syms x y z phi theta psi
disp('Матрица g(n):')
disp(vpa(g([x; y; z; phi; theta; psi]), 4))

%% РАСЧЕТ ДИНАМИКИ МАНИПУЛЯТОРА %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ЗАГРУЗИТЬ ПАРАМЕТРЫ ИЗ ФАЙЛА
manip_param

% Homogeneous transformations solution
Ti = cell(n,1);

A = @(a, alpha, d, theta) ...
        ([cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha) a*cos(theta); ...
        sin(theta) cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta); ...
        0 sin(alpha) cos(alpha) d; ...
        0 0 0 1]);

T0 = eye(4); T0(1:3,1:3) = R0; T0(1:3,4) = origin;
T = T0;
for i = 1:n
    temp = A(DH(i,1),DH(i,2),DH(i,3),DH(i,4));
    T = T*temp;
    Ti{i} = simplify(T);
end

%% Центры масс
com = cell(n,1);
for i=1:n    
P = Ti{i}*[[1;0;0;0] [0;1;0;0] [0;0;1;0] [c{i};1]];
    x = P(1,4);
    y = P(2,4);
    z = P(3,4);
    com{i} = [x; y; z];  
end

%% Inertia matrix and kinetic energy
qd = sym('qd', [n 1], 'real'); % "q dot" - the first derivative of the q's in time (joint velocities)
 
% Velocity Jacobians
Jv = cell(1,n);
Jw = cell(1,n);
for i = 1:n
    z = T0(1:3,3);
    o = T0(1:3,4);
    Jvt = sym(zeros(3,n));
    Jwt = sym(zeros(3,n));
    for j = 1:i
        Jvt(:,j) = cross(z, com{i}-o);
        Jwt(:,j) = z;
        z = Ti{j}(1:3,3);
        o = Ti{j}(1:3,4);
    end
    Jv{i} = Jvt;
    Jw{i} = Jwt;
end

%% ВТОРОЙ СПОСОБ (РЕЗУЛЬТАТ ТОТ ЖЕ!) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jv = cell(1,n);
% Jw = cell(1,n);
% for i = 1:n
%     z = T0(1:3,3);
%     Jvt = sym(zeros(3,n));
%     Jwt = sym(zeros(3,n));
%     for j = 1:i
%         Jvt(:,j) = diff(com{i}, q(j));
%         Jwt(:,j) = z;
%         z = Ti{j}(1:3,3);
%     end
%     Jv{i} = Jvt;
%     Jw{i} = Jwt;
% end

%% Potential energy
% Potential energy solution
P = eye(4);
PE = 0;
% M_sym = Inertia matrix solution & PE = Poterntial Energy
M_sym = 0;
for i = 1:n
    R = Ti{i}(1:3,1:3);
    M_sym = M_sym + (m(i)*Jv{i}'*Jv{i} + Jw{i}'*R*I{i}*R'*Jw{i});
    PE = PE + m(i)*g_accel*com{i}(3);
end

%% Equations of motion
qdd = sym('qdd', [n 1], 'real'); % "q double dot" - the second derivative of the q's in time (joint accelerations)

% The Christoffel symbols
c = zeros(n,n,n,'sym');
for k = 1:n
    for i = 1:n
        for j = 1:n
            c(i,j,k) = 0.5 * (diff(M_sym(k,j),q(i)) + diff(M_sym(k,i),q(j)) - diff(M_sym(i,j),q(k)));
        end
    end
end

% The coriolis matrix
C_sym = zeros(n,n,'sym');
for k = 1:n
    for j = 1:n
        temp = 0;
        for i = 1:n
            temp = temp + c(i,j,k)*qd(i);
        end
        C_sym(k,j) = temp; 
    end
end

% The gravitation terms
g_sym = zeros(n,1,'sym');
for k = 1:n
    g_sym(k) = diff(PE,q(k));
end

% M_sym = simplify(M_sym);
% C_sym = simplify(C_sym);
% g_sym = simplify(g_sym);

matlabFunction(M_sym,'File','get_M','Vars',{[q]});
matlabFunction(C_sym,'File','get_C','Vars',{[q;qd]});
matlabFunction(g_sym,'File','get_g','Vars',{[q]});

%% МОДЕЛИРОВАНИЕ ДВИЖЕНИЯ
eta0 = [0, 0, -5, 0, 0, 0]'; % начальная глубина 5 метров
v0 = [0, 0, 0, 0, 0, 0]';
q0 = deg2rad([0, 0]'); % начальные углы
dq0 = [0, 0]';

tau = @(t)(heaviside(t).*[10 0 0 0 0 0 0 0]' - heaviside(t-15).*[10 0 0 0 0 0 0 0]'  );
t_end = 30; dt = 0.01;
[t,Y] = ode45(@(t,y)odefcn(t,y,M,C,D,g,tau), 0:dt:t_end, [eta0; v0; q0; dq0]);

%% ПОСТРОЕНИЯ ГРАФИКОВ
% Графики для AUV
figure(1), clf
v = Y(:,7:12);
subplot(2,2,1), title('Скорости (линейные)'), hold on, grid on
plot(t, v(:,1:3)), xlabel('t, сек'), ylabel('Скорость, м/c'), xlim([0 t_end])
legend('u(t)', 'v(t)', 'w(t)', 'Location', 'Best')
subplot(2,2,2), title('Скорости (угловые)'), hold on, grid on
plot(t, v(:,4:end)), xlabel('t, сек'), ylabel('Скорость, рад/c'), xlim([0 t_end])
legend('p(t)', 'q(t)', 'r(t)', 'Location', 'Best')
eta = Y(:,1:6);
subplot(2,2,3), title('Положения (по осям)'), hold on, grid on
plot(t, eta(:,1:3)), xlabel('t, сек'), ylabel('Положения, м'), xlim([0 t_end])
legend('x(t)', 'y(t)', 'z(t)', 'Location', 'Best')
subplot(2,2,4), title('Положения (углы Эйлера)'), hold on, grid on
plot(t, eta(:,4:end)), xlabel('t, сек'), ylabel('Положения, рад'), xlim([0 t_end])
legend('\phi(t)', '\theta(t)', '\psi(t)', 'Location', 'Best')

figure(2), clf, plot3(eta(:,1), eta(:,2), eta(:,3)), title('Позиционирование AUV')
hold on
plot3(eta(end,1), eta(end,2), eta(end,3), 'rO') 
legend('траектория AUV', 'конечная точка', 'Location', 'Best')
xlabel('x(t)'), ylabel('y(t)'), zlabel('z(t)'), grid on
set(gca, 'YDir', 'reverse');
view([-1,1,1])
axis equal

figure(3), clf, plot3(eta(:,1), eta(:,2), eta(:,3)), title('Позиционирование AUV')
hold on
plot3(eta(end,1), eta(end,2), eta(end,3), 'rO') 
xlabel('x(t)'), ylabel('y(t)'), zlabel('z(t)'), grid on
set(gca, 'YDir', 'reverse');
view([-1,1,1])
axis equal
legend('off')

rectangular_plot(eta(1,:), L*10^-3, W*10^-3, H*10^-3, '--r')
rectangular_plot(eta(end,:), L*10^-3, W*10^-3, H*10^-3, 'r')

q_vect = Y(:,13:end-n);
hold on
manip_plot(eta(1,:)', q, q_vect(1,:)', com, Ti, l, r, '--r')
hold on
manip_plot(eta(end,:)', q, q_vect(end,:)', com, Ti, l, r, 'r')


% Графики для манипулятора
figure(4), clf
title('Значения q_i(t)'), hold on, grid on
plot(t, rad2deg(q_vect)), xlabel('t, сек'), ylabel('Положение, deg'), xlim([0 t_end])
labels = [];
for i = 1:n
    labels = [labels; strcat('q_', num2str(i), '(t)')];
end
legend(cellstr(labels))