function [Fc,lambda_trial] = AugmentedContactForce(ContactBody,TargetBody,approach)

    ContactPointfunc = approach.ContactPointfunc;
    penalty = approach.penalty;
    lambda_old = approach.lambda.meaning;

    Fcont = zeros(ContactBody.nx,1);
    Ftarg = zeros(TargetBody.nx,1);

    [ContactPoints,ContactPointsElements,ContactAreas] =ContactPointfunc(ContactBody);

    numberOfPoints = size(ContactPoints,1);

    if isempty(lambda_old)
        lambda_old = zeros(numberOfPoints,1);
    end

    lambda_trial = zeros(numberOfPoints,1);

    for i = 1:numberOfPoints

        ContactPoint = ContactPoints(i,:); 

        % We need the signed gap for every candidate point.
        Outcome = FindPoint(TargetBody,ContactPoint);

        % Projection onto the target surface is unavailable.
        if ~isfield(Outcome,"Index")
            lambda_trial(i) = 0;
            continue
        end

        signedGap = Outcome.Gap;
        lambda_trial(i) = max(0, lambda_old(i) - penalty*ContactAreas(i)*signedGap);

        if lambda_trial(i) == 0
            continue
        end

        Normal = Outcome.Normal;

        % Existing force-direction convention.
        Normal_cont = -Normal;
        Normal_targ =  Normal;

        % Contact element.
        element_cont = ContactPointsElements(i);
        [X_cont,U_cont] = GetCoorDisp(element_cont,ContactBody.nloc,ContactBody.P0,ContactBody.u);
        [xi_cont,eta_cont] = FindIsoCoord(X_cont,U_cont,ContactPoint');

        % Target element.
        element_targ = Outcome.Index;
        [X_targ,U_targ] = GetCoorDisp(element_targ,TargetBody.nloc,TargetBody.P0,TargetBody.u);
        [xi_targ,eta_targ] = FindIsoCoord(X_targ,U_targ,Outcome.Position);

        Fcont_loc =  lambda_trial(i)*Normal_cont;
        Ftarg_loc =  lambda_trial(i)*Normal_targ;

        DOFpositions_cont = ContactBody.xloc(element_cont,:);
        DOFpositions_targ = TargetBody.xloc(element_targ,:);

        Fcont(DOFpositions_cont) = Fcont(DOFpositions_cont)+ Nm_2412(xi_cont,eta_cont)'*Fcont_loc;
        Ftarg(DOFpositions_targ) = Ftarg(DOFpositions_targ)+ Nm_2412(xi_targ,eta_targ)'*Ftarg_loc;
    end

    Fc = [Fcont;Ftarg];
end