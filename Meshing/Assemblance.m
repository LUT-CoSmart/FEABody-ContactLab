function [Body1, Body2, uu_bc, deltaf, lambda] = Assemblance(iteration,Solution,Body1,Body2,DofsFunction,Stiffness,approach,lambdaBackTrack)
    
    persistent Bn_m1 ff_bc_old uu_bc_old 

    Type = approach.Type;
    Name = approach.Name;
    lambda = approach.lambda.meaning; % keeping just for standard Lagrange method
    penalty = approach.penalty;
    approach.lambda.converge = true;
    bc = [Body1.bc Body2.bc]; % total logical vector of constrains
    Fext = [Body1.Fext.vec; Body2.Fext.vec]; % Assemblance of external forces
           
    Fe = [Body1.Fint.vec; Body2.Fint.vec]; % Assemblance of elastic forces    
    ff =  Fe - Fext; % residuals
    ff_bc = ff(bc);
    Fext_bc = Fext(bc);
    
    % Assemblance of elastic stiffness
    Ke = [            Body1.Fint.K zeros(Body1.nx,Body2.nx);
          zeros(Body2.nx,Body1.nx)            Body2.Fint.K];
    % Removal of the fixed dofs
    Ke_bc = Ke(bc,bc); 
        
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
            ff_bc = [ff_bc; zeros(m,1)]; % this works for the cases where initial gap is zero  
            
        otherwise
            ff_bc = ff_bc + DofsFunction(bc);            
            K_bc = Ke_bc + Stiffness(bc,bc);
            
    end    
    
    if Solution == "Newton-Rapson" || (Solution == "Newton-Broyden" && iteration == 1)
        num = cond(K_bc);
        Bn_m1 = inv(K_bc);
        if num > 1e12
           D = diag(1./sqrt(sum(K_bc.^2,2)));          
           Bn_m1 = (D*K_bc)\D;
        else 
           Bn_m1 = inv(K_bc);           
        end
         uu_bc = - Bn_m1*ff_bc;
    elseif Solution == "Newton-Broyden" 
       zn = ff_bc - ff_bc_old;
       s = uu_bc_old;
       denom = s' * Bn_m1 * zn;
       Bn_m1 = Bn_m1 + (s - Bn_m1 * zn ) * s' * Bn_m1 ./ denom ; 
       uu_bc = - Bn_m1*ff_bc;

    else
        error('Unkown solutiuon method');
    end

    % for Newton-Broyden:
    ff_bc_old = ff_bc;
    uu_bc_old = lambdaBackTrack * uu_bc;
    
    deltaf = ff_bc(1:size(Ke_bc,1))/norm(Fext_bc);
    uu_bc = lambdaBackTrack *uu_bc(1:Body1.ndof + Body2.ndof); % removing potential additional DOFs from Lagrange-based approaches
    
    % Displacement separation
    Body1.u(Body1.bc) = Body1.u(Body1.bc) + uu_bc(1:Body1.ndof);
    Body2.u(Body2.bc) = Body2.u(Body2.bc) + uu_bc(Body1.ndof + 1:end);

   