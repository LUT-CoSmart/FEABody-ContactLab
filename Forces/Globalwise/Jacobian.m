function K= Jac_stepwise(Body1,Body2,AimFunction)

    nx = Body1.nx + Body2.nx;
    y0 = AimFunction(Body1,Body2);
    K = zeros(length(y0),nx);  
    I_vec=zeros(nx,1); 
    h = 2*sqrt(eps);
    % Backup original coordinates
    u1_backup = Body1.u;
    u2_backup = Body2.u;
    for ii = 1:nx
        % start variationing 
        I_vec(ii)=1;
        % this split is to distribute coord. between bodies 
        Body1.u = u1_backup + h*I_vec(1:Body1.nx); 
        Body2.u = u2_backup + h*I_vec(1+Body1.nx:end);
        y_p = AimFunction(Body1,Body2);
        K(:,ii) = (y_p - y0) / h;
        I_vec(ii)=0;  
    end    
    % returning all back to normal
    Body1.u = u1_backup;
    Body2.u = u2_backup;

