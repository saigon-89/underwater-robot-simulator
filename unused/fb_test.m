%% МОДЕЛИРОВАНИЕ ДВИЖЕНИЯ
eta0 = [0, 0, -15, 0, 0, 0]'; % начальная глубина 15 метров
v0 = [0, 0, 0, 0, 0, 0]'; 

u = [10, 20, -20, 0, 0, pi/4]'; % уставка

t_end = 1500; dt = 0.01;

global forces
forces = [];
[t,Y] = ode45(@(t,y)odefcn_fb(t,y,M,C,D,g,u), 0:dt:t_end, [eta0; v0]);
figure
for i=1:6
    plot(linspace(0,t_end,length(forces)), forces(:,i)), grid on
    hold on
end
legend('\tau_1(t)','\tau_2(t)','\tau_3(t)','\tau_4(t)','\tau_5(t)','\tau_6(t)')
xlabel('t, сек'), ylabel('Силы и моменты')
xlim([0 t_end])

%% ПОСТРОЕНИЯ ГРАФИКОВ
figure
v = Y(:,7:end);
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

figure, plot3(eta(:,1), eta(:,2), eta(:,3)), title('Позиционирование AUV')
hold on
plot3(eta(end,1), eta(end,2), eta(end,3), 'gO') 
plot3(u(1), u(2), u(3), 'rO') 
legend('траектория AUV', 'конечная точка', 'уставка', 'Location', 'Best')
xlabel('x(t)'), ylabel('y(t)'), zlabel('z(t)'), grid on
set(gca, 'YDir', 'reverse');
view([-1,1,1])
axis equal

legend('off')
rectangular_plot(eta(1,:), L*10^-3, W*10^-3, H*10^-3, '--r')
rectangular_plot(eta(end,:), L*10^-3, W*10^-3, H*10^-3, 'r')
