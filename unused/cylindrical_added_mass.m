function A = cylindrical_added_mass(m, r, L, rho)

    if nargin < 4
       rho = 1000; 
    end
    
    A = zeros(6,6);

    A(1,1) = 0.1*m;
    A(2,2) = pi*rho*r^2*L*10^-3;
    A(3,3) = A(2,2);
    A(4,4) = 0;
    A(5,5) = pi*rho*r^2*L/12;
    A(6,6) = A(5,5);
    
end

