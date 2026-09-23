function [DofsFunction,Stiffness] = Contact(Body1,Body2,approach,step,Solution)

    DOFsNumber = Body1.nx + Body2.nx;
    DofsFunction = zeros(DOFsNumber,1);
    Stiffness = zeros(numel(DofsFunction),DOFsNumber);

    if approach.Type ~= "None"
        
        AimFunction = @(Body1,Body2) approach.AimFunction(Body1,Body2,approach);
        DofsFunction = AimFunction(Body1,Body2);

        if step == 1 || Solution == "Newton-Rapson" 
            
            switch approach.perturbation   
                case "automatic"
                    DofsFunction_y = @(t,y) AimFunctionWrapper(t,y,Body1,Body2,AimFunction);
                    Stiffness = numjac(DofsFunction_y,0,[Body1.u(:);Body2.u(:)], DofsFunction,1e-3,[]);
    
                case "incremental"
                    Stiffness = Jacobian(Body1,Body2,AimFunction);
            end
        end    
    end