function K = Jac_stepwise_Im(Body1,Body2,AimFunction)

    nx = Body1.nx + Body2.nx;

    % Evaluate once to determine the output size
    y0 = AimFunction(Body1,Body2);
    K = zeros(length(y0),nx);

    % Complex-step size
    h = 1e-30;

    % Backup original displacement vectors
    u1_backup = Body1.u;
    u2_backup = Body2.u;

    for ii = 1:nx

        % Start from the original real displacement vectors
        Body1.u = complex(u1_backup);
        Body2.u = complex(u2_backup);

        % Add an imaginary perturbation to one DOF
        if ii <= Body1.nx
            Body1.u(ii) = Body1.u(ii) + 1i*h;
        else
            jj = ii - Body1.nx;
            Body2.u(jj) = Body2.u(jj) + 1i*h;
        end

        % Evaluate the residual with the perturbed displacement
        y_trial = AimFunction(Body1,Body2);

        % Corresponding Jacobian column
        K(:,ii) = imag(y_trial)/h;
    end

    % Restore the original displacement vectors
    Body1.u = u1_backup;
    Body2.u = u2_backup;

end