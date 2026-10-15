function Sigma = Sigma_raw_2412_Neo(E,nu,U,X,xi,eta)

    mu = E/(2*(1+nu));
    lambda   = E*nu/((1+nu)*(1-2*nu));

    F = F_2412(U,X,xi,eta);
    J = det(F);
    if real(J) <= 0
        error('Sigma_raw_2412_Neo: element inverted (J = %g)', real(J));
    end
    C_inv = inv(F.'*F);

    Sigma = mu*(eye(2) - C_inv) + lambda*log(J)*C_inv;
end