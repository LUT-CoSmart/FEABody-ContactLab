function [Fc,lambda_trial] = AugmentedNitscheContactForce(ContactBody,TargetBody,approach)
    ContactPointfunc = approach.ContactPointfunc;
    Fcont = zeros(ContactBody.nx,1);
    Ftarg = zeros(TargetBody.nx,1);
    
    [ContactPoints,ContactPointsElements] = ContactPointfunc(ContactBody);
    numberOfPoints = size(ContactPoints,1);
    area = ContactBody.Lx/ContactBody.nElems.x*ContactBody.Lz/approach.numberOfPoints;
    gamma = area/approach.penalty;

    lambda_old = approach.lambda.meaning;
    if isempty(lambda_old) || length(lambda_old) ~= numberOfPoints
        lambda_old = zeros(numberOfPoints,1);
    end
    lambda_trial = lambda_old;

    for i = 1:numberOfPoints
        ContactPoint = ContactPoints(i,:);
        Outcome = FindPoint(TargetBody,ContactPoint);

        if ~isfield(Outcome,"Index")
            lambda_trial(i) = 0;
            continue
        end

        signedGap = Outcome.Gap; % it is negative for penetration
        Normal = Outcome.Normal(:);
        Normal_cont = -Normal;
        Normal_targ = Normal;
       
        % contact
        element_cont = ContactPointsElements(i);
        [X_cont,U_cont] = GetCoorDisp(element_cont,ContactBody.nloc,ContactBody.P0,ContactBody.u);
        [xi_cont,eta_cont] = FindIsoCoord(X_cont,U_cont,ContactPoint');
        Nm_cont = Nm_2412(xi_cont,eta_cont);
        Sigma_cont = Sigma_2412(ContactBody.E,ContactBody.nu,U_cont,X_cont,xi_cont,eta_cont);
        sigma_nn_cont = Normal_cont.'*Sigma_cont*Normal_cont;
        dsigma_cont = NormalStressDerivative_Im(ContactBody.E,ContactBody.nu,U_cont,X_cont,xi_cont,eta_cont,Normal_cont);
        
        % target
        element_targ = Outcome.Index;
        [X_targ,U_targ] = GetCoorDisp(element_targ,TargetBody.nloc,TargetBody.P0,TargetBody.u);
        [xi_targ,eta_targ] = FindIsoCoord(X_targ,U_targ,Outcome.Position);
        Nm_targ = Nm_2412(xi_targ,eta_targ);
        Sigma_targ = Sigma_2412(TargetBody.E,TargetBody.nu,U_targ,X_targ,xi_targ,eta_targ);
        sigma_nn_targ = Normal_targ.'*Sigma_targ*Normal_targ;
        dsigma_targ = NormalStressDerivative_Im(TargetBody.E,TargetBody.nu,U_targ,X_targ,xi_targ,eta_targ,Normal_targ);

        sigma_nn = 0.5*(sigma_nn_cont+sigma_nn_targ); % likely it will be positive
        dsigma_cont = 0.5*dsigma_cont;
        dsigma_targ = 0.5*dsigma_targ;
        
        % Positive compressive pressure; negative gap means penetration.
        augmentedPressure = max(0,lambda_old(i)-signedGap/gamma);
        
        % Relax only the outer multiplier update.
        lambda_trial(i) = augmentedPressure;
        
        % sigma_nn is tension-positive:
        % pressure-stress agreement means augmentedPressure + sigma_nn = 0.
        stressMismatch = augmentedPressure + sigma_nn;
        
        Fcont_loc = augmentedPressure*(Nm_cont.'*Normal_cont) - gamma*stressMismatch*dsigma_cont;       
        Ftarg_loc = augmentedPressure*(Nm_targ.'*Normal_targ) - gamma*stressMismatch*dsigma_targ;
        
        % Redistribution over the nodes
        DOFpositions_cont = ContactBody.xloc(element_cont,:);
        DOFpositions_targ = TargetBody.xloc(element_targ,:);
        
        Fcont(DOFpositions_cont) = Fcont(DOFpositions_cont)+Fcont_loc*area;% *area here improves convergence
        Ftarg(DOFpositions_targ) = Ftarg(DOFpositions_targ)+Ftarg_loc*area;


    end

    Fc = [Fcont;Ftarg];

end