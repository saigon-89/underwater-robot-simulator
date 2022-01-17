function dy = odefcn_fb(t,y,M,C,D,g,u)
    global forces

    % Matrix form
    % A = [ zeros(6), eye(6);
    %       zeros(6), - M\C(y(7:end)) - M\D(y(7:end)) ];
    % B = [ zeros(6); 
    %       M\eye(6) ];
    % dy = A*y + B*err - [ zeros(6,1); M\g(y(1:6)) ];

    % Controller
    K = [10 1 1 1 1 1]';
    
    % Calculate error and control signal
    err = (u - y(1:6));
    ctrl = K .* err;
    
    % Saturation
    % hlim = [10 10 10 1 1 1]; 
    % llim = -hlim;
    % for i=1:length(err)
    %    if(err(i) > hlim(i)) 
    %        err(i) = hlim(i); 
    %    else if(err(i) < llim(i))
    %        err(i) = llim
    %    end
    % end
    
    forces = [forces; err'];
    
    dy = zeros(12,1);
    dy(1:6) = y(7:end);
    dy(7:end) = - M\C(y(7:end))*y(7:end) - M\D(y(7:end))*y(7:end) - ...
       M\g(y(1:6)) + M\ctrl;

end