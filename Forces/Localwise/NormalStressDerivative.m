function dsigma = NormalStressDerivative(E,nu,U,X,xi,eta,Normal)

    h = 2*sqrt(eps);
    

    % Normal stress at the original displacement vector
    Sigma_0 = Sigma_2412(E,nu,U,X,xi,eta);   
    sigmaNN_0 = Normal.'*Sigma_0*Normal;

    dsigma = zeros(8,1);
    I_vec = zeros(8,1);

    for j = 1:8
        % Perturb one displacement component
        I_vec(j) = 1;
        Uh = U - h*I_vec;

        Sigmah = Sigma_2412(E,nu,Uh,X,xi,eta);        
        sigmaNN_h = Normal.'*Sigmah*Normal;
        dsigma(j) = (sigmaNN_0 - sigmaNN_h)/h;

        I_vec(j) = 0;
    end

end