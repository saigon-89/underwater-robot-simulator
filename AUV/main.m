%% GUI ÔÎÐÌÀ ÄËß ÂÂÎÄÀ ÏÀÐÀÌÅÒÐÎÂ AUV 
prompt = {'Äëèíà ROV [mm]:',
    'Âûñîòà ROV [mm]:',
    'Øèðèíà ROV [mm]:',
    'Ïëîòíîñòü æèäêîñòè [kg/m^3]:',  
    'Ïëîùàäü ôðîíòîâîé ïðîåêöèè [mm^2]:', 
    'Ïëîùàäü áîêîâîé ïðîåêöèè [mm^2]:',
    'Ïëîùàäü âåðõíåé ïðîåêöèè [mm^2]:',
    'Ìàññà ROV [kg]:',
    'Âåêòîð r^g_c [mm]: ',
    'Âåêòîð r^b_c [mm]: '};
dlgtitle = 'Âõîäíûå ïàðàìåòðû';
dims = [1 45];
definput = {'0','0','0','1000','0','0','0','0','[0, 0, 0]', '[0, 0, 0]'};
options.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,options);

L = str2num(answer{1});
H = str2num(answer{2});
W = str2num(answer{3});
rho = str2num(answer{4});
PF = str2num(answer{5});
PS = str2num(answer{6});
PT = str2num(answer{7});
m = str2num(answer{8});
r_g_c = str2num(answer{9}); % âåêòîð îò íà÷àëà êîîðäèíàò äî öåíòðà òÿæåñòè
r_b_c = str2num(answer{10}); % âåêòîð îò íà÷àëà êîîðäèíàò äî ãåîìåòðè÷åñêîãî öåíòðà

%% ÕÀÐÀÊÒÅÐÍÎ ÄËß ÏÐßÌÎÓÃÎËÜÍÎÉ ÏÐÈÇÌÛ
I0  = diag(m*[(W^2+H^2) (L^2+H^2) (W^2+L^2)]/12).*10^-6; % òåíçîð èíåðöèè (1.6)

%% ÂÑÏÎÌÎÃÀÒÅËÜÍÛÅ ÔÓÍÊÖÈÈ
% ïðåîáðàçîâàíèå â êîñîñèììåòðè÷íóþ ìàòðèöó
S = @(x)([ 0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0 ]); % (1.5)
r_g_c = r_g_c / 10^3; % ïåðåâîäèì [mm] â [m]
r_b_c = r_b_c / 10^3; % ïåðåâîäèì [mm] â [m]

%% ÐÀÑ×ÅÒ ÌÀÒÐÈÖÛ M
% ðàñ÷åò M_RB
M_RB = [ m*eye(3) -m*S(r_b_c); m*S(r_b_c) I0 ]; % (1.4)
% ðàñ÷åò M_A
M_A = rectangular_added_mass(L, H, W, rho, PF, PS, PT); % [6]
% ðàñ÷åò M
M = M_RB + M_A; % (1.3)

%% ÐÀÑ×ÅÒ C(v)
% ðàñ÷åò C_RB(v)
C_RB = @(v)([ zeros(3) -m*S(v(1:3)); -m*S(v(1:3)) -S(diag(I0).*v(4:end)) ]); % (1.7)
% ðàñ÷åò C_A(v)
diag_M_A = diag(M_A);
C_A = @(v)([ zeros(3) -S(diag_M_A(1:3).*v(1:3)); % (1.8)
    -S(diag_M_A(1:3).*v(1:3)) -S(diag_M_A(4:end).*v(4:end))]);
% ðàñ÷åò C(v)
C = @(v)(C_RB(v) + C_A(v)); % (1.3)

%% ÐÀÑ×ÅÒ g(n)
B = m * 9.81; % (1.9)
g = @(eta)([ 0; 
    0; 
    0;
    -(r_g_c(2)*B-r_b_c(2)*B)*cos(eta(5))*cos(eta(4)) + (r_g_c(3)*B-r_b_c(3)*B)*cos(eta(5))*sin(eta(4));
    (r_g_c(3)*B-r_b_c(3)*B)*sin(eta(5)) + (r_g_c(1)*B-r_b_c(1)*B)*cos(eta(5))*cos(eta(4));
    -(r_g_c(1)*B-r_b_c(1)*B)*cos(eta(5))*sin(eta(4)) - (r_g_c(2)*B-r_b_c(2)*B)*sin(eta(5)) ]); % (1.10)

%% ÐÀÑ×ÅÒ D(v)
D = @(v)(rectangular_damping(L, H, W, rho, PF, PS, PT, M_RB, M_A, B, r_g_c, r_b_c).*v); % [6]

%% ÂÛÂÎÄ ÐÅÇÓËÜÒÀÒÎÂ
disp('Ìàòðèöà M:')
disp(vpa(M, 4))

syms u v w p q r
disp('Ìàòðèöà Ñ(v):')
disp(vpa(C([u; v; w; p; q; r]), 4))

disp('Ìàòðèöà D(v):')
disp(vpa(D([u; v; w; p; q; r]), 4))

syms x y z phi theta psi
disp('Ìàòðèöà g(n):')
disp(vpa(g([x; y; z; phi; theta; psi]), 4))

%% ÌÎÄÅËÈÐÎÂÀÍÈÅ ÄÂÈÆÅÍÈß
eta0 = [0, 0, -5, 0, 0, 0]'; % íà÷àëüíàÿ ãëóáèíà 5 ìåòðîâ
v0 = [0, 0, 0, 0, 0, 0]';
tau = @(t)(heaviside(t).*ones(6,1));
t_end = 120; dt = 0.001;
[t,Y] = ode45(@(t,y)odefcn(t,y,M,C,D,g,tau), 0:dt:t_end, [eta0; v0]);

%% ÏÎÑÒÐÎÅÍÈß ÃÐÀÔÈÊÎÂ
figure
v = Y(:,7:end);
subplot(2,2,1), title('Ñêîðîñòè (ëèíåéíûå)'), hold on, grid on
plot(t, v(:,1:3)), xlabel('t, ñåê'), ylabel('Ñêîðîñòü, ì/c'), xlim([0 t_end])
legend('u(t)', 'v(t)', 'w(t)', 'Location', 'Best')
subplot(2,2,2), title('Ñêîðîñòè (óãëîâûå)'), hold on, grid on
plot(t, v(:,4:end)), xlabel('t, ñåê'), ylabel('Ñêîðîñòü, ðàä/c'), xlim([0 t_end])
legend('p(t)', 'q(t)', 'r(t)', 'Location', 'Best')
eta = Y(:,1:6);
subplot(2,2,3), title('Ïîëîæåíèÿ (ïî îñÿì)'), hold on, grid on
plot(t, eta(:,1:3)), xlabel('t, ñåê'), ylabel('Ïîëîæåíèÿ, ì'), xlim([0 t_end])
legend('x(t)', 'y(t)', 'z(t)', 'Location', 'Best')
subplot(2,2,4), title('Ïîëîæåíèÿ (óãëû Ýéëåðà)'), hold on, grid on
plot(t, eta(:,4:end)), xlabel('t, ñåê'), ylabel('Ïîëîæåíèÿ, ðàä'), xlim([0 t_end])
legend('\phi(t)', '\theta(t)', '\psi(t)', 'Location', 'Best')

figure, plot3(eta(:,1), eta(:,2), eta(:,3)), title('Ïîçèöèîíèðîâàíèå AUV')
hold on
plot3(eta(end,1), eta(end,2), eta(end,3), 'rO') 
legend('òðàåêòîðèÿ AUV', 'êîíå÷íàÿ òî÷êà', 'Location', 'Best')
xlabel('x(t)'), ylabel('y(t)'), zlabel('z(t)'), grid on
set(gca, 'YDir', 'reverse');
view([-1,1,1])
axis equal
