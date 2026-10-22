function [ContactPointPositions, ContactTraction] = ContactPressureAtPoints(ContactBody, TargetBody, approach)
% Contact pressure t_N of the chosen method at every contact point (zero where there is no contact).
% Same formulas as the contact-force functions; post-processing only (not for augmented Nitsche).
    [ContactPointPositions, ContactPointsElements, ContactAreas] = approach.ContactPointfunc(ContactBody);
    numberOfPoints = size(ContactPointPositions,1);
    ContactTraction = zeros(numberOfPoints,1);

    PointwisePenalty = approach.penalty*ones(numberOfPoints,1);      % penalty-Nitsche keeps its own p_n per point
    if startsWith(approach.Name, "penalty-Nitsche") && numel(approach.lambda.meaning) == numberOfPoints
        RaisedPoints = approach.lambda.meaning > 0;
        PointwisePenalty(RaisedPoints) = approach.lambda.meaning(RaisedPoints);
    end

    for i = 1:numberOfPoints
        Outcome = FindPoint(TargetBody, ContactPointPositions(i,:));
        if ~isfield(Outcome,"Index")
            continue
        end
        penetration = -Outcome.Gap;                                   % Gap is negative for penetration
        switch approach.Name
            case "Penalty"
                ContactTraction(i) = approach.penalty*max(0, penetration);
            case "Augumented Lagrange"
                ContactTraction(i) = max(0, approach.lambda.meaning(i) + approach.penalty*ContactAreas(i)*penetration)/ContactAreas(i);
            otherwise                                                  % Nitsche, penalty-Nitsche: stress-based pressure
                Normal = Outcome.Normal(:);
                element_cont = ContactPointsElements(i);
                [X_cont,U_cont] = GetCoorDisp(element_cont,ContactBody.nloc,ContactBody.P0,ContactBody.u);
                [xi_cont,eta_cont] = FindIsoCoord(X_cont,U_cont,ContactPointPositions(i,:)');
                Sigma_cont = Sigma_2412(ContactBody.E,ContactBody.nu,U_cont,X_cont,xi_cont,eta_cont,ContactBody.StressFunction);
                element_targ = Outcome.Index;
                [X_targ,U_targ] = GetCoorDisp(element_targ,TargetBody.nloc,TargetBody.P0,TargetBody.u);
                [xi_targ,eta_targ] = FindIsoCoord(X_targ,U_targ,Outcome.Position);
                Sigma_targ = Sigma_2412(TargetBody.E,TargetBody.nu,U_targ,X_targ,xi_targ,eta_targ,TargetBody.StressFunction);
                sigma_nn = 0.5*(Normal.'*Sigma_cont*Normal + Normal.'*Sigma_targ*Normal);   % n'*sigma*n is the same for n and -n
                ContactTraction(i) = max(0, -sigma_nn + PointwisePenalty(i)*penetration);
        end
    end
end