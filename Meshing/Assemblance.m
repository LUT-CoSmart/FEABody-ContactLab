function [Body1, Body2, uu_bc, deltaf, lambda] = Assemblance(Solution,Body1,Body2,DofsFunction,Stiffness,approach)
    
    Type = approach.Type;
    Name = approach.Name;
    lambda = approach.lambda.meaning; % keeping just for standard Lagrange method
    penalty = approach.penalty;
    approach.lambda.converge = true;
    bc = [Body1.bc Body2.bc]; % total logical vector of constrains
    Fext = [Body1.Fext.vec; Body2.Fext.vec]; % Assemblance of external forces
    
    % Assemblance of elastic stiffness
    Ke = [            Body1.Fint.K zeros(Body1.nx,Body2.nx);
          zeros(Body2.nx,Body1.nx)            Body2.Fint.K];
    
    Fe = [Body1.Fint.vec; Body2.Fint.vec]; % Assemblance of elastic forces
    
    ff =  Fe - Fext; % residuals

    % Removal of the fixed dofs
    Ke_bc = Ke(bc,bc);  
    ff_bc = ff(bc);
    Fext_bc = Fext(bc);
    
    % Methods
    switch Type
        case "Lagrange"
            m = size(Stiffness,1);
            Stiffness_bc = Stiffness(:,bc);

            switch Name 
                case "perturbed Lagrange"
                    Matrix_add = -1/penalty * ones(m);

                case "Lagrange" 
                    Matrix_add = zeros(m);
            end

            K_bc = [        Ke_bc Stiffness_bc';
                     Stiffness_bc  Matrix_add]; 
            ff_bc = [ff_bc; zeros(m,1)];
            
        otherwise % "Penalty" and "None"
            ff_bc = ff_bc + DofsFunction(bc);
            K_bc = Ke_bc + Stiffness(bc,bc);
            
            if Name == "Augumented Lagrange"
                ff_bc = ff_bc + lambda; % contributions from the updated contact forces after the previous iteration
                lambda = lambda + DofsFunction(bc); 
            end

    end    
    
    if Solution == "Newton-Rapson"
        num = cond(K_bc);
        if num > 1e12
           D = diag(1./sqrt(sum(K_bc.^2,2)));
           uu_bc =- (D*K_bc)\(D*ff_bc);
        else 
            uu_bc = -K_bc\ff_bc;
        end
    else
        error('Unkown solutiuon method');
    end
    deltaf = ff_bc(1:size(Ke_bc,1))/norm(Fext_bc);
    uu_bc = uu_bc(1:Body1.ndof + Body2.ndof); % removing potential additional DOFs from Lagrange-based approaches
    
    % Displacement separation
    Body1.u(Body1.bc) = Body1.u(Body1.bc) + uu_bc(1:Body1.ndof);
    Body2.u(Body2.bc) = Body2.u(Body2.bc) + uu_bc(Body1.ndof + 1:end);

   