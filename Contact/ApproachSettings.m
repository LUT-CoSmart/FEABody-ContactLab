function [approach,backtrack] = ApproachSettings(approachBasis,approachSubtype,ContactPointfunc,...
                                GapfuncPairs,Perturbation,PointsofInterest,ContactBody,TargetBody,alpha)
    
    backtrack.meaning = false;% staring back track for the best solution to find an equlibrium
    approach.lambda.meaning = [];
    approach.lambda.imax = 1; % activation maximal iterations for augumemted algorithms 
    approach.perturbation = Perturbation;
    approach.ContactPointfunc = ContactPointfunc;
    approach.numberOfPoints = PointsofInterest.n;
    

    if approachSubtype == "Augumented Lagrange" || approachSubtype == "penalty-Nitsche" || approachSubtype == "Augumented Nitsche"  
        [InitialContactPoints,~] = ContactPointfunc(ContactBody);
        approach.lambda.meaning = zeros(size(InitialContactPoints,1),1);
        approach.lambda.imax = 100;
    end
    
    % h  = min(ContactBody.Lx/ContactBody.nElems.x, TargetBody.Lx/TargetBody.nElems.x);
    hC = min(vecnorm(diff(ContactBody.P0(ContactBody.contact.nodalid,:)),2,2));
    hT = min(vecnorm(diff(TargetBody.P0(TargetBody.contact.nodalid,:)),2,2));
    h  = min(hC,hT); % propsed by Claude - do compare all for any orientation 
    
    approach.penalty = alpha*max(ContactBody.E,TargetBody.E)/h;  

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
            AimFunction = @(Body1,Body2,approach) AugmentedContactForce(Body1,Body2,approach);

        elseif approachSubtype == "penalty-Nitsche"
            backtrack.meaning = true;
            approach.gapTolerance = 1e-4*h;
            AimFunction = @(Body1,Body2,approach) PenaltyNitscheContactForce(Body1,Body2,approach);
            
        elseif approachSubtype == "Augumented Nitsche"          
            AimFunction = @(Body1,Body2,approach) AugmentedNitscheContactForce(Body1,Body2,approach);

        elseif approachSubtype == "Penalty"
            backtrack.meaning = true;
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
    
    % Activate backtracking algorithm
    if backtrack.meaning 
       backtrack.lambdaList  = [1.0, 0.5, 0.25, 0.125];   
    else
       backtrack.lambdaList = 1;
    end
   
