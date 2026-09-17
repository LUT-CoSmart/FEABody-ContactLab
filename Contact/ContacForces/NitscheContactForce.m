function Fc = NitscheContactForce(ContactBody,TargetBody,approach)
    ContactPointfunc = approach.ContactPointfunc;
    Fcont = zeros(ContactBody.nx,1);
    Ftarg = zeros(TargetBody.nx,1);
    
    [ContactPoints,ContactPointsElements] = ContactPointfunc(ContactBody);
    numberOfPoints = size(ContactPoints,1);
    area = ContactBody.Lx/ContactBody.nElems.x*ContactBody.Lz/approach.numberOfPoints;
    gamma = area/approach.penalty;

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

        element_cont = ContactPointsElements(i);
        [X_cont,U_cont] = GetCoorDisp(element_cont,ContactBody.nloc,ContactBody.P0,ContactBody.u);
        [xi_cont,eta_cont] = FindIsoCoord(X_cont,U_cont,ContactPoint');
        Nm_cont = Nm_2412(xi_cont,eta_cont);
        Sigma_cont = Sigma_2412(ContactBody.E,ContactBody.nu,U_cont,X_cont,xi_cont,eta_cont);
        sigma_nn_cont = Normal_cont.'*Sigma_cont*Normal_cont;
        dsigma_cont = NormalStressDerivative(ContactBody.E,ContactBody.nu,U_cont,X_cont,xi_cont,eta_cont,Normal_cont);
        
        element_targ = Outcome.Index;
        [X_targ,U_targ] = GetCoorDisp(element_targ,TargetBody.nloc,TargetBody.P0,TargetBody.u);
        [xi_targ,eta_targ] = FindIsoCoord(X_targ,U_targ,Outcome.Position);
        Nm_targ = Nm_2412(xi_targ,eta_targ);
        Sigma_targ = Sigma_2412(TargetBody.E,TargetBody.nu,U_targ,X_targ,xi_targ,eta_targ);
        sigma_nn_targ = Normal_targ.'*Sigma_targ*Normal_targ;
        dsigma_targ = NormalStressDerivative(TargetBody.E,TargetBody.nu,U_targ,X_targ,xi_targ,eta_targ,Normal_targ);
        
        sigma_nn = 0.5*(sigma_nn_cont+sigma_nn_targ); % likely it will be positive
        dsigma_cont = 0.5*dsigma_cont;
        dsigma_targ = 0.5*dsigma_targ;
        
        penetration = -signedGap;
        dpenetration_cont = Nm_cont.'*Normal_cont;
        dpenetration_targ = Nm_targ.'*Normal_targ;
        Pgamma = penetration-gamma*sigma_nn;
        PgammaPositive = max(0,Pgamma);
        
        Rcont = -gamma*sigma_nn*dsigma_cont+(PgammaPositive/gamma)*(dpenetration_cont-gamma*dsigma_cont);
        Rtarg = -gamma*sigma_nn*dsigma_targ+(PgammaPositive/gamma)*(dpenetration_targ-gamma*dsigma_targ);
        
        DOFpositions_cont = ContactBody.xloc(element_cont,:);
        DOFpositions_targ = TargetBody.xloc(element_targ,:);
        
        Fcont(DOFpositions_cont) = Fcont(DOFpositions_cont)+Rcont*area;
        Ftarg(DOFpositions_targ) = Ftarg(DOFpositions_targ)+Rtarg*area;

    end

    Fc = [Fcont;Ftarg];
end
