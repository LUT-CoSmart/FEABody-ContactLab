function approach = ApproachSettings(approachBasis,approachSubtype,ContactPointfunc, GapfuncPairs, Perturbation,backtrack)
    
    approach.lambda.meaning = [];
    approach.penalty =1e9;  
    approach.perturbation = Perturbation;
    approach.ContactPointfunc = ContactPointfunc; % for augumented lagrange
    
    

    % sanity check, that approach and subtype are correlating 
    if approachBasis == "Lagrange"
       backtrack = false;
       warning("backtrack is off for this set up");
       allowed = ["Lagrange","perturbed Lagrange"];  

       if ~any(approachSubtype == allowed)
          warning("Invalid approachSubtype for approachBasis='Lagrange, substituted to Lagrange'");
          approachSubtype = "Lagrange";
       end

       approach.Name = approachSubtype; 
       AimFunction = GapfuncPairs; 

    elseif approachBasis == "Penalty"
        
        allowed = ["Penalty", "Nitshe-linear", "Nitshe-nonlinear", "Nitshe-nonlinear-all", "Augumented Lagrange"];
        
             
        if ~any(approachSubtype == allowed)
          warning("Invalid approachSubtype for approachBasis='Penalty, substituted to Penalty");
          approachSubtype = "Penalty";          
        end
        
        approach.Name = approachSubtype; 

        if approachSubtype == "Augumented Lagrange"
            backtrack = false;
            warning("backtrack is off for this set up");
            approach.penalty = 1e7;  % decreasing parameter for better stability, method operates with any small penalty 
            AimFunction = @(Body1,Body2) AugmentedContactForce(Body1,Body2,approach);
        else
            AimFunction = @(Body1,Body2) ContactForce(Body1,Body2,approach);
        end
                    
    else
        warning("Contact is not activated");
        AimFunction = [];
        approachBasis = "None"; 
        approach.Name = "None";
    end  

    approach.Type = approachBasis;   
    approach.AimFunction = AimFunction;
    
    if backtrack 
       approach.lambdaList  = [1.0, 0.5, 0.25, 0.125];   
    else
       approach.lambdaList = 1;
    end