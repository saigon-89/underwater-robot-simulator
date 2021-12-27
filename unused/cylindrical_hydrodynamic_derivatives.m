function vect = cylindrical_hydrodynamic_derivatives(m, r, L, rho)

    if nargin < 4
       rho = 1000; 
    end
    
    X_du = -0.1*m;
    Y_dv = -pi*rho*r^2*L;
    Z_dw = Y_dv;
    K_dp = 0;
    M_dq = -pi*rho*r^2*L^-3;
    N_dr = M_dq;
    vect = [X_du; Y_dv; Z_dw; K_dp; M_dq; N_dr];
    
end

