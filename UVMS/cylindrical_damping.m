function [BL, BQ] = cylindrical_damping(L, D, rho, PF, PS, PT, M_RB, M_A, B, r_g_c, r_b_c)
    
    %% INPUT VALUES
    lambda = 0.16; % scaling linear/quadratic
    I44 = M_RB(4,4); % Moment of inertia in roll
    I55 = M_RB(5,5); % Moment of inertia in pitch
    A44 = M_A(4,4); % Added mass in roll
    A55 = M_A(5,5); % Added mass in pitch
    C = abs(B*(r_g_c(3)-r_b_c(3))); % Restoring coefficient (pitch=roll)
    
    %% COEFFICENTS
    CpXY = PT/(L*D*pi/2); % Projected Area Coefficient XY
    CpYZ = PF/(pi*D^2/4); % Projected Area Coefficient YZ
    CpZX = PS/(L*D*pi/2); % Projected Area Coefficient XZ
    
    %% DRAG COEFFICIENTS (3D)
    Data3D_par = [0,1.15;0.5,1.10;1.0,0.93;1.5,0.85;2.0,0.83;3.0,0.85;4.0,0.85;5.0,0.85];
    Drag3D_par=spline(Data3D_par(:,1),Data3D_par(:,2));
    Data3D_vert = [1.0,0.64;1.98,0.68;2.98,0.74;10.0,0.82;20.0,0.91;40.0,0.98;999,1.2];
    Drag3D_vert=spline(Data3D_vert(:,1),Data3D_vert(:,2));

    %% NONLINEAR DAMPING
    % Final Surge nonlinear damping
    LD=L/D;
    BQ(1,1)=0.5*rho*ppval(Drag3D_par,(LD))*(pi*D^2/4)*10^-9*CpYZ;
    % Sway
    LD=L/D;
    BQ(2,2)=rho*0.5*ppval(Drag3D_vert,(LD))*(L*D*pi/2)*10^-6*CpZX;
    % Heave
    LD=L/D;
    BQ(3,3)=rho*0.5*ppval(Drag3D_vert,(LD))*(L*D*pi/2)*10^-6*CpXY;
    
    %% ROLL
    LD=L/D;
    Fh=rho*(1/6)*ppval(Drag2D,(LD))*(H/2)*L*10^-6*CpZX;
    Mh=Fh*(3/4)*((H/2)*10^-3)^3;
    LD=L/D;
    Fv=rho*(1/6)*ppval(Drag2D,(LD))*(W/2)*L*10^-6*CpXY;
    Mv=Fv*(3/4)*((W/2)*10^-3)^3;
    BQ(4,4)=(2*Mv+2*Mh);
    
    %% PITCH
    LD=L/D;
    Fh=rho*(1/6)*ppval(Drag2D,(LD))*(H/2)*W*10^-6*CpYZ;
    Mh=Fh*(3/4)*((H/2)*10^-3)^3;
    LD=L/D;
    Fv=rho*(1/6)*ppval(Drag2D,(LD))*(L/2)*W*10^-6*CpXY;
    Mv=Fv*(3/4)*((L/2)*10^-3)^3;
    BQ(5,5)=(2*Mv+2*Mh);

    %% YAW
    LD=L/D;
    Fh=rho*(1/6)*ppval(Drag2D,(LD))*(W/2)*H*10^-6*CpYZ;
    Mh=Fh*(3/4)*((W/2)*10^-3)^3;
    LD=L/D;
    Fv=rho*(1/6)*ppval(Drag2D,(LD))*(L/2)*H*10^-6*CpZX;
    Mv=Fv*(3/4)*((L/2)*10^-3)^3;
    BQ(6,6)=(2*Mv+2*Mh);
    
    %% LINEAR VISCOUS DAMPING
    % Roll and Pitch
    BL(4,4)=2*0.025*(I44+A44)*sqrt(C/(I44+A44));
    BL(5,5)=2*0.025*(I55+A55)*sqrt(C/(I55+A55));
    lambda0=0.16;
    lambda1=BL(5,5)/BQ(5,5);
    % Surge, sway, heave and yaw
    BL(1,1)=BQ(1,1)*lambda0;
    BL(2,2)=BQ(2,2)*lambda0;
    BL(3,3)=BQ(3,3)*lambda0;
    BL(6,6)=BQ(6,6)*lambda1;
    
end

