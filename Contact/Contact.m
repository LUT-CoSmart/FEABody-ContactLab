function [DofsFunction,Stiffness] = Contact(Body1,Body2,approach,step,Solution)

    DOFsNumber = Body1.nx + Body2.nx;
    DofsFunction = zeros(DOFsNumber,1);
    Stiffness = zeros(DOFsNumber);

    if approach.Type ~= "None"

        %% TODO: combine these two together    
        if approach.Name == "Augumented Lagrange"             
            % lambda_old from Current input of approach
            AimFunction = @(Body1,Body2) AugmentedContactForce(Body1,Body2,approach); 

        elseif approach.Name == "penalty-Nitsche"             
            % lambda_old from Current input of approach
            AimFunction = @(Body1,Body2) PenaltyNitscheContactForce(Body1,Body2,approach);     

        elseif approach.Name == "Augmented Nitsche"
            AimFunction = @(Body1,Body2) AugmentedNitscheContactForce(Body1,Body2,approach);    
        
        else
            AimFunction = approach.AimFunction;
        end

        DofsFunction = AimFunction(Body1,Body2);

        if step == 1 || Solution == "Newton-Rapson" 
            
            switch approach.perturbation   
                case "automatic"
                    DofsFunction_y = @(t,y) AimFunctionWrapper(t,y,Body1,Body2,AimFunction);
                    Stiffness = numjac(DofsFunction_y,0,[Body1.u(:);Body2.u(:)], DofsFunction,1e-3,[]);
    
                case "incremental"
                    Stiffness = Jac_stepwise(Body1,Body2,AimFunction);
            end
        end    
    end