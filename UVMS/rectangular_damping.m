function [BL, BQ] = rectangular_damping(L, H, W, rho, PF, PS, PT, M_RB, M_A, B, r_g_c, r_b_c)
    
    %% INPUT VALUES
    lambda = 0.16; % scaling linear/quadratic
    I44 = M_RB(4,4); % Moment of inertia in roll
    I55 = M_RB(5,5); % Moment of inertia in pitch
    A44 = M_A(4,4); % Added mass in roll
    A55 = M_A(5,5); % Added mass in pitch
    C = abs(B*(r_g_c(3)-r_b_c(3))); % Restoring coefficient (pitch=roll)
    
    %% COEFFICENTS
    CpXY = PT/(L*W); % Projected Area Coefficient XY
    CpYZ = PF/(H*W); % Projected Area Coefficient YZ
    CpZX = PS/(L*H); % Projected Area Coefficient XZ
    
    %% DRAG COEFFICIENTS (2D)
    Data2D = [0.5,2.5;1.5,1.8;2.5,1.4;6,0.89]; 
    Drag2D = spline(Data2D(:,1),Data2D(:,2));

    %% DRAG COEFFICIENTS (3D)
    Data3D =[0,1.25;0.5,1.25;1,1.15;1.5,0.97;2,0.87;...
        2.5,0.9;3,0.93;4,0.95;5,0.95];
    Drag3D=spline(Data3D(:,1),Data3D(:,2));

    %% NONLINEAR DAMPING
    % Surge 3D
    LD=L/((H+W)/2);
    BQ3D=ppval(Drag3D,(LD));
    % Surge 2D
    LD=L/W;
    BQ2D=ppval(Drag2D,(LD));
    
    lambda = BQ3D/BQ2D;
    % Final Surge nonlinear damping
    LD=L/((H+W)/2);
    BQ(1,1)=0.5*rho*ppval(Drag2D,(LD))*H*W*10^-6*CpYZ*lambda;
    % Sway
    LD=W/H;
    BQ(2,2)=rho*0.5*ppval(Drag2D,(LD))*L*H*10^-6*CpZX*lambda;
    % Heave
    LD=H/W;
    ppval(Drag2D,(LD));
    BQ(3,3)=rho*0.5*ppval(Drag2D,(LD))*L*W*10^-6*CpXY*lambda;
    
    %% ROLL
    LD=W/(H/2);
    Fh=rho*(1/6)*ppval(Drag2D,(LD))*(H/2)*L*10^-6*CpZX*lambda;
    Mh=Fh*(3/4)*((H/2)*10^-3)^3;
    LD=H/(W/2);
    Fv=rho*(1/6)*ppval(Drag2D,(LD))*(W/2)*L*10^-6*CpXY*lambda;
    Mv=Fv*(3/4)*((W/2)*10^-3)^3;
    BQ(4,4)= (2*Mv+2*Mh);
    
    %% PITCH
    LD=L/(H/2);
    Fh=rho*(1/6)*ppval(Drag2D,(LD))*(H/2)*W*10^-6*CpYZ*lambda;
    Mh=Fh*(3/4)*((H/2)*10^-3)^3;
    LD=H/(L/2);
    Fv=rho*(1/6)*ppval(Drag2D,(LD))*(L/2)*W*10^-6*CpXY*lambda;
    Mv=Fv*(3/4)*((L/2)*10^-3)^3;
    BQ(5,5)=(2*Mv+2*Mh);

    %% YAW
    LD=L/(W/2);
    Fh=rho*(1/6)*ppval(Drag2D,(LD))*(W/2)*H*10^-6*CpYZ*lambda;
    Mh=Fh*(3/4)*((W/2)*10^-3)^3;
    LD=W/(L/2);
    Fv=rho*(1/6)*ppval(Drag2D,(LD))*(L/2)*H*10^-6*CpZX*lambda;
    Mv=Fv*(3/4)*((L/2)*10^-3)^3;
    BQ(6,6)= (2*Mv+2*Mh);
    
    %% LINEAR VISCOUS DAMPING
    % Roll and Pitch
    BL(4,4)= 2*0.025*(I44+A44)*sqrt(C/(I44+A44));
    BL(5,5)=2*0.025*(I55+A55)*sqrt(C/(I55+A55));
    lambda0=0.16;
    lambda1=BL(5,5)/BQ(5,5);
    % Surge, sway, heave and yaw
    BL(1,1)=BQ(1,1)*lambda0;
    BL(2,2)=BQ(2,2)*lambda0;
    BL(3,3)=BQ(3,3)*lambda0;
    BL(6,6)=BQ(6,6)*lambda1;
    
end

