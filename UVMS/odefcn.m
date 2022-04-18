function dy = odefcn(t,y,M,C,D,g,tau)
% M - AUV inertia matrix
% C - AUV Coriolis matrix
% D - AUV damping matrix
% g - AUV gravity vector
% tau - UVMS forces/moments vector
% 
% Y is the state vector = [eta v];
% v is 6-vector of velocities:
% - u(t), v(t), w(t) - Linear Velocities, m/sec
% - p(t), q(t), r(t) - Angular Velocities, rad/sec
% 
% eta is 6-vector of Positions:
% - x(t), y(t), z(t) - Linear Positions, m
% - \phi(t), \theta(t), \psi(t) - Positions (Euler angles), rad
% 
% [t,Y] = ode45(@(t,y)odefcn(t,y,M,C,D,g,tau), 0:dt:t_end, [eta0; v0]);
% v = Y(:,7:end);
% eta = Y(:,1:6);
%

eta = y(1:6); 
v = y(7:12); 
q = y(13:14); 
dq = y(15:16); 
n = (numel(y) - 12) / 2; % manipulator links

%     M_full = [ M,          zeros(6,n);
%                zeros(n,6), get_M([q]) ];
%     C_full = [ C(v),       zeros(6,n);
%                zeros(n,6), get_C([q; dq]) ];
%     D_full = [ D(v),       zeros(6,n);
%                zeros(n,6), 0.5*eye(n) ];
%     g_full = [ g(eta);
%                get_g([q]) ];

temp = tau(t);
tau_auv = temp(1:6);
tau_manip = temp(7:end);

dy = zeros(16,1);
dy(1:6) = y(7:12);
dy(7:12) = M\( -C(v)*v - D(v)*v - g(eta) + tau_auv);
dy(13:14) = y(15:16);
dy(15:16) = get_M([q]) \ ( -get_C([q; dq])*dq - 0.9*eye(2)*dq - get_g(q) + tau_manip - [tau_auv(1);0]);
%     dy = zeros(12+2*n,1);
%     dy(1:6+n) = y(6+n+1:end);
%     dy(6+n+1:end) = M_full \ ( -C_full*[v;dq] - D_full*[v;dq] - g_full + tau(t) );


    
%% STATE-SPACE TERMS   
%     A = [ zeros(6), eye(6);
%           zeros(6), M\( -C_full*[v;dq] - D_full*[v;dq] - g_full ) ];
%     B = [ zeros(6); 
%           M\eye(6) ];
%     dy = A*y + B*tau(t);
end
