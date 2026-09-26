function Sigma = Sigma_2412(E,nu,U,X,xi,eta) % stress recover,
    % such name because the whole program was written for raw stress function

    [gp2,~] = gauleg2(-1,1,2);
    gp2 = sort(gp2);
    [gp1,~] = gauleg2(-1,1,1);
         
    gm = gp2(1);
    gp = gp2(2);
    
    % Gauss point 1
    S1 = Sigma_raw_2412(E,nu,U,X,gm,gm);
    F1 = F_2412(U,X,gm,gm);
    C1 = 1/det(F1) * F1 * S1 * F1.'; % . in vcase of the usage of Im_ based jacobian
    
    % Gauss point 2
    S2 = Sigma_raw_2412(E,nu,U,X,gp,gm);
    F2 = F_2412(U,X,gp,gm);
    C2 = 1/det(F2) * F2 * S2 * F2.';
    
    % Gauss point 3
    S3 = Sigma_raw_2412(E,nu,U,X,gp,gp);
    F3 = F_2412(U,X,gp,gp);
    C3 = 1/det(F3) * F3 * S3 * F3.';
    
    % Gauss point 4
    S4 = Sigma_raw_2412(E,nu,U,X,gm,gp);
    F4 = F_2412(U,X,gm,gp);
    C4 = 1/det(F4) * F4 * S4 * F4.';

    xi_r  = xi*sqrt(3); % standart scale factor for Q4
    eta_r = eta*sqrt(3);

    Nm = Nm_2412(xi_r,eta_r);

    N1 = Nm(1,1);
    N2 = Nm(1,3);
    N3 = Nm(1,5);
    N4 = Nm(1,7);

    Sigma = zeros(2);

    Sigma(1,1) = N1*C1(1,1) + N2*C2(1,1) + N3*C3(1,1) + N4*C4(1,1);
    Sigma(2,2) = N1*C1(2,2) + N2*C2(2,2) + N3*C3(2,2) + N4*C4(2,2);

    Scenter = Sigma_raw_2412(E,nu,U,X,gp1,gp1);
    Fcenter = F_2412(U,X,gp1,gp1);

    Ccenter = 1/det(Fcenter) * Fcenter * Scenter * Fcenter.';

    Sigma(1,2) = Ccenter(1,2);
    Sigma(2,1) = Ccenter(2,1);

end