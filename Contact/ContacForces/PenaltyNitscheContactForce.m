function [Fc,pn_trial,maxPenetration] = PenaltyNitscheContactForce(ContactBody,TargetBody,approach)
    ContactPointfunc = approach.ContactPointfunc;
    Fcont = zeros(ContactBody.nx,1);
    Ftarg = zeros(TargetBody.nx,1);

    [ContactPoints,ContactPointsElements,ContactAreas] = ContactPointfunc(ContactBody); 
    numberOfPoints = size(ContactPoints,1); 
    maxPenetration = 0; 
    gapTolerance = approach.gapTolerance; 

    pn_old = approach.lambda.meaning;
    if isempty(pn_old) || numel(pn_old) ~= numberOfPoints
        pn_old = approach.penalty*ones(numberOfPoints,1);
    else
        pn_old(pn_old <= 0) = approach.penalty;
    end
    pn_trial = pn_old;

    for i = 1:numberOfPoints 
        ContactPoint = ContactPoints(i,:); 
        Outcome = FindPoint(TargetBody,ContactPoint); 

        if ~isfield(Outcome,"Index") 
            continue 
        end 

        signedGap = Outcome.Gap;        
        Normal = Outcome.Normal(:); 
        Normal_cont = -Normal; 
        Normal_targ =  Normal;
       
         % contact
        element_cont = ContactPointsElements(i);
        [X_cont,U_cont] = GetCoorDisp(element_cont,ContactBody.nloc,ContactBody.P0,ContactBody.u);
        [xi_cont,eta_cont] = FindIsoCoord(X_cont,U_cont,ContactPoint');
        Nm_cont = Nm_2412(xi_cont,eta_cont);
        Sigma_cont = Sigma_2412(ContactBody.E,ContactBody.nu,U_cont,X_cont,xi_cont,eta_cont);
        sigma_nn_cont = Normal_cont.'*Sigma_cont*Normal_cont;
        dsigma_cont = NormalStressDerivative(ContactBody.E,ContactBody.nu,U_cont,X_cont,xi_cont,eta_cont,Normal_cont);
        
        % target
        element_targ = Outcome.Index;
        [X_targ,U_targ] = GetCoorDisp(element_targ,TargetBody.nloc,TargetBody.P0,TargetBody.u);
        [xi_targ,eta_targ] = FindIsoCoord(X_targ,U_targ,Outcome.Position);
        Nm_targ = Nm_2412(xi_targ,eta_targ);
        Sigma_targ = Sigma_2412(TargetBody.E,TargetBody.nu,U_targ,X_targ,xi_targ,eta_targ);
        sigma_nn_targ = Normal_targ.'*Sigma_targ*Normal_targ;
        dsigma_targ = NormalStressDerivative(TargetBody.E,TargetBody.nu,U_targ,X_targ,xi_targ,eta_targ,Normal_targ);

        
        sigma_nn = 0.5*(sigma_nn_cont + sigma_nn_targ);
        dsigma_cont = 0.5*dsigma_cont;
        dsigma_targ = 0.5*dsigma_targ;
        
        penetration = max(0, -signedGap);
        
        maxPenetration = max(maxPenetration, penetration);
        

        if penetration > gapTolerance
            pn_target_gap = pn_old(i)*penetration/gapTolerance;
            pn_target_sigma = max(0, -sigma_nn)/gapTolerance;        
            pn_target = max(pn_target_gap, pn_target_sigma);        
            pn_trial(i) = min(10*pn_old(i), max(pn_old(i), pn_target));
        end
                
        gamma = 1/pn_old(i);
               
        pressure = max(0, -sigma_nn - signedGap/gamma);
                
        Fcont_loc = ContactAreas(i)*(pressure*(Nm_cont.'*Normal_cont) - gamma*(sigma_nn + pressure)*dsigma_cont);        
        Ftarg_loc = ContactAreas(i)*(pressure*(Nm_targ.'*Normal_targ) - gamma*(sigma_nn + pressure)*dsigma_targ);

        DOFpositions_cont = ContactBody.xloc(element_cont,:); 
        DOFpositions_targ = TargetBody.xloc(element_targ,:); 
       
        Fcont(DOFpositions_cont) = Fcont(DOFpositions_cont)+Fcont_loc;% *area here improves convergence
        Ftarg(DOFpositions_targ) = Ftarg(DOFpositions_targ)+Ftarg_loc;

    end 

    Fc = [Fcont;Ftarg]; 
