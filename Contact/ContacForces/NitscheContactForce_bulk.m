function [Fc,pn_trial] = NitscheContactForce_bulk(ContactBody,TargetBody,approach)

    ContactPointfunc = approach.ContactPointfunc; 
    pn_old = approach.lambda.meaning; 
    gapTolerance = approach.gapTolerance; 
    Fcont = zeros(ContactBody.nx,1); 
    Ftarg = zeros(TargetBody.nx,1); 

    
    [ContactPoints,ContactPointsElements] = ContactPointfunc(ContactBody); 
    numberOfPoints = size(ContactPoints,1); 
    
    if isempty(pn_old) || length(pn_old) ~= numberOfPoints 
        pn_old = ones(numberOfPoints,1)*approach.penalty; 
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
        Normal_targ = Normal;

        area = ContactBody.Lx/ContactBody.nElems.x*ContactBody.Lz;

        % Contact element 
        element_cont = ContactPointsElements(i); 
        [X_cont,U_cont] = GetCoorDisp(element_cont,ContactBody.nloc,ContactBody.P0,ContactBody.u); 
        [xi_cont,eta_cont] = FindIsoCoord(X_cont,U_cont,ContactPoint'); 
        Sigma_cont = Sigma_2412(ContactBody.E,ContactBody.nu,U_cont,X_cont,xi_cont,eta_cont);
        sigma_nn_cont = Normal_cont.'*Sigma_cont*Normal_cont; 
        xivec_cont = [xi_cont;eta_cont]; 
        nabla_Sigma_nn_cont = NablaSigma_NN(ContactBody.E,ContactBody.nu,U_cont,X_cont,xivec_cont,Normal_cont);
        nabla_sigma_cont = 0.5*nabla_Sigma_nn_cont;
        
        % Target element 
        element_targ = Outcome.Index; 
        [X_targ,U_targ] = GetCoorDisp(element_targ,TargetBody.nloc,TargetBody.P0,TargetBody.u); 
        [xi_targ,eta_targ] = FindIsoCoord(X_targ,U_targ,Outcome.Position); 
        Sigma_targ = Sigma_2412(TargetBody.E,TargetBody.nu,U_targ,X_targ,xi_targ,eta_targ); 
        sigma_nn_targ = Normal_targ.'*Sigma_targ*Normal_targ; xivec_targ = [xi_targ;eta_targ];
        nabla_Sigma_nn_targ = NablaSigma_NN(TargetBody.E,TargetBody.nu,U_targ,X_targ,xivec_targ,Normal_targ); 
        nabla_sigma_targ = 0.5*nabla_Sigma_nn_targ; 
        
        % Integrated compression-positive bulk stress 
        sigma = 0.5*(sigma_nn_cont+sigma_nn_targ)*area; 

        % Target penalty and current Nitsche contribution 
        if sigma > 0 
            pn_target = abs(sigma/gapTolerance); % we assume that we hace compression and sigma is "-" 
        else % just in case 
            pn_target = pn_old(i); 
        end 

        % Projected contact force 
        pressure = max(0,-pn_old(i)*signedGap); % as soon as we reach desired condtion sigma-pn*gap is fullfiled automatically  from the condition of the target gap 
        if pressure == 0 
            continue
        end 
        
        % Penalty proposed for the next outer iteration 
        if pn_old(i) < pn_target 
            pn_trial(i) = min(10*pn_old(i),pn_target); 
        end 
        
        % Local Nitsche force vectors 
         
        Fcont_loc = pressure*Normal_cont + signedGap*nabla_sigma_cont*area; 
        Ftarg_loc = pressure*Normal_targ + signedGap*nabla_sigma_targ*area; 
                      
        DOFpositions_cont = ContactBody.xloc(element_cont,:); 
        DOFpositions_targ = TargetBody.xloc(element_targ,:); 
       
        Fcont(DOFpositions_cont) = Fcont(DOFpositions_cont)+Nm_2412(xi_cont,eta_cont)'*Fcont_loc; 
        Ftarg(DOFpositions_targ) = Ftarg(DOFpositions_targ)+Nm_2412(xi_targ,eta_targ)'*Ftarg_loc; 

    end 

    Fc = [Fcont;Ftarg]; 
end
