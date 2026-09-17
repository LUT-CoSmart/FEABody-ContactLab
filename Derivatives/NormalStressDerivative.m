function dsigma = NormalStressDerivative(E,nu,U,X,xi,eta,Normal)

    % test the idea of imagery things
    %% TODO: make it automatic as well
    h = 1e-30;
    dsigma = zeros(length(U),1);
    for j = 1:length(U)
        U_trial = complex(U);
        U_trial(j) = U_trial(j)+1i*h;
        Sigma_trial = Sigma_2412(E,nu,U_trial,X,xi,eta);
        dsigma(j) = imag(Normal.'*Sigma_trial*Normal)/h;
    end
end