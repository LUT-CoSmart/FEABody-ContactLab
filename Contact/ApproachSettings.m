function [approach,backtrack] = ApproachSettings(approachBasis,approachSubtype,ContactPointfunc,GapfuncPairs,Perturbation,PointsofInterest)
    
    backtrack.meaning = false;% staring back track for the best solution to find an equlibrium
    approach.lambda.meaning = [];
    approach.lambda.imax = 10;
    % approach.penalty =1e9;  
    approach.perturbation = Perturbation;
    approach.ContactPointfunc = ContactPointfunc; % for augumented lagrange
    approach.numberOfPoints = PointsofInterest.n;

    % sanity check, that approach and subtype are correlating 
    if approachBasis == "Lagrange"
       allowed = ["Lagrange","perturbed Lagrange"];  

       if ~any(approachSubtype == allowed)
          warning("Invalid approachSubtype for approachBasis='Lagrange, substituted to Lagrange'");
          approachSubtype = "Lagrange";
       end

       approach.Name = approachSubtype; 
       AimFunction = GapfuncPairs; 

    elseif approachBasis == "Penalty"
        
        allowed = ["Penalty", "Nitsche", "penalty-Nitsche", "Augumented Nitsche", "Augumented Lagrange"];
                             
        if ~any(approachSubtype == allowed)
          warning("Invalid approachSubtype for approachBasis='Penalty, substituted to Penalty");
          approachSubtype = "Penalty";          
        end
        
        approach.Name = approachSubtype; 

        if approachSubtype == "Augumented Lagrange"
            approach.penalty = 1e7;  % decreasing parameter for better stability, method operates with any small penalty 
            AimFunction = @(Body1,Body2,approach) AugmentedContactForce(Body1,Body2,approach);

        elseif approachSubtype == "penalty-Nitsche"
            backtrack.meaning = true;
            approach.gapTolerance = 2e-4;
            AimFunction = @(Body1,Body2,approach) PenaltyNitscheContactForce(Body1,Body2,approach);
            
        elseif approachSubtype == "Augumented Nitsche"          
            AimFunction = @(Body1,Body2,approach) AugmentedNitscheContactForce(Body1,Body2,approach);

        elseif approachSubtype == "Penalty"
            % backtrack.meaning = true;
            backtrack.meaning = false;
            AimFunction = @(Body1,Body2,approach) PenaltyContactForce(Body1,Body2,approach);
        
        elseif approachSubtype == "Nitsche"
            backtrack.meaning = true;
            AimFunction = @(Body1,Body2,approach) NitscheContactForce(Body1,Body2,approach);
        end
                    
    else
        warning("Contact is not activated");
        AimFunction = [];
        approachBasis = "None"; 
        approach.Name = "None";
    end  

    approach.Type = approachBasis;   
    approach.AimFunction = AimFunction;
    
    if backtrack.meaning 
       backtrack.lambdaList  = [1.0, 0.5, 0.25, 0.125];   
    else
       backtrack.lambdaList = 1;
    end