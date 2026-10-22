%clc, clear, close all;
if ~exist('SWEEP','var')|| ~SWEEP % we have apper file to deal with
    clc, clear, close all
    SWEEP = false;
else
    clear Body1 Body2 approach backtrack Result Assemblance BestStep Contact
end

format long;
addpath(genpath(pwd));

%########## Bodies Data : Element positioning (from (0.0) coord. ) ########
% Body 1
Body1.Lx = 0.01;
Body1.Ly = 0.01;
Body1.Lz = 0.1;
Body1.E = 2e11;
Body1.nu = 0.3;
Radius = 0.05;
Body1 = SetMaterial(Body1, "KS");

% Body 2
Body2.Lx = 0.04;
Body2.Ly = 0.02;
Body2.Lz = 0.1;
Body2.E = 2e11;
Body2.nu = 0.3;
Body2 = SetMaterial(Body2, "KS");


Body1.shift.x = (Body2.Lx-Body1.Lx)/2;
Body1.shift.y = Body2.Ly - 1e-6;        % the lowest point of the indenter touches the block at x = 0

%#################### Mesh #########################################
dx1 = 20;
dy1 = 8;

dx2 = 40;
dy2 = 10;

Body1.nElems.x = dx1;
Body1.nElems.y = dy1;

Body2.nElems.x = dx2;
Body2.nElems.y = dy2;

Body1 = CreateFEMesh(Body1);
Body2 = CreateFEMesh(Body2);

%#################### BC ###########################################
Body1.loc.x = 'all';
Body1.loc.y = Body1.Ly;
Body1 = CreateBC(Body1,'all');

Body2.loc.x = 'all';
Body2.loc.y = 0;
Body2 = CreateBC(Body2,'all');

%##################### Prescribed displacement ##############
VerticalDisplacement = 2.5e-4;
AppliedDisp.x = 'all';
AppliedDisp.y = Body1.Ly;
AppliedDispNodes = FindGlobNodalID(Body1.P0, AppliedDisp, Body1.shift);
AppliedDispVertical = xlocChosen(Body1.DofsAtNode, AppliedDispNodes, 2);

%##################### Contact surfaces #############################
Body1.contact.loc.x = 'all';
Body1.contact.loc.y = 0;
Body1.contact.nodalid = FindGlobNodalID(Body1.P0,Body1.contact.loc,Body1.shift);

Body2.contact.loc.x = 'all';
Body2.contact.loc.y = Body2.Ly;
Body2.contact.nodalid = FindGlobNodalID(Body2.P0,Body2.contact.loc,Body2.shift);

%##################### Making bottom surface round #############################
Body1 = RoundingBottom(Body1, Radius);

% shape check
% Visualization(Body1);
% Visualization(Body2);

%##################### Contact ######################################
% Subtypes
% Penalty: Penalty, Augumented Lagrange,
%          Nitsche, penalty-Nitsche,  Augumented Nitsche
% Lagrange: Lagrange, perturbed Lagrange
if ~SWEEP
    approachBasis = "Penalty";
    approachSubtype = "penalty-Nitsche";
    alpha = 100;
    %##################### Newton iter. parameters ######################
    imax = 50;
    tol = 1e-4;
    LoadType = "linear"; % no external forces here; the indentation is applied linearly in the steps
    steps = 1;
    Solution = "Newton-Rapson"; % Options: Newton-Rapson, Newton-Broyden
    PenetrationTol = 1e-3;
    ShowVisualization = true;
    ShowNodeNumbers = false;
    WhatoToShow = "sigma_VM"; % options: "ux", "uy", "u_total", "sigma_xx", "sigma_yy", "sigma_xy", "sigma_VM"
end

PointsofInterest.Name = "Gauss";
PointsofInterest.n = 2;         

[ContactPointfunc, Gapfunc, GapfuncPairs]  = ContactPointSetting(PointsofInterest,approachBasis);
Perturbation = "automatic"; % Options: "automatic", "incremental"

Backtrack = false;
[approach,backtrack] = ApproachSettings(approachBasis,approachSubtype,ContactPointfunc,...
                       GapfuncPairs,Perturbation,PointsofInterest,Body1,Body2,alpha,PenetrationTol,Backtrack);
if ~SWEEP
    approach.theta = -1; % assymetric nitsche
else
    approach.theta = theta;
end

%#################### Processing ######################
total_steps = 0;
titertot=0;
Niter = 0;
failed = false;

for ii = 1:steps

        approach.lambda.step = 0;
        approach.lambda.converge = false;

        Body1 = CreateFext(ii,steps,Body1,LoadType);  % these give zeros, but necessary
        Body2 = CreateFext(ii,steps,Body2,LoadType);  % these give zeros, but necessary
        Body1.u(AppliedDispVertical) = -VerticalDisplacement*ii/steps; % prescribed disp

        while (~approach.lambda.converge) && (approach.lambda.step < approach.lambda.imax)% special case for Augumented Lagrange
                approach.lambda.step = approach.lambda.step + 1;
                for lambdaBackTrack = backtrack.lambdaList % backtrack loop

                    if lambdaBackTrack < 1
                        fprintf("!!! Starting  Backtrack for lambda = %f !!!! \n",lambdaBackTrack)
                    end

                    for jj = 1:imax % contact convergence
                        tic;
                        Niter = Niter + 1;
                        % interaction of two bodies
                        [DofsFunction, Stiffness] = Contact(Body1,Body2,approach,jj,Solution);

                        Gap = Gapfunc(Body1,Body2);
                        % inner forces of the each body
                        Body1 = Elastic(Body1);
                        Body2 = Elastic(Body2);

                        [Body1, Body2, uu_bc, deltaf, lambda_next] = Assemblance(jj,Solution,Body1, Body2, DofsFunction,Stiffness,approach,lambdaBackTrack);

                        titer=toc;
                        titertot=titertot+titer;

                        conv = printStatus(approachBasis, deltaf, uu_bc, tol, ii, jj, imax, steps, titertot, Gap);
                        if conv
                            break;
                        end

                        [Body1,Body2] = BestStep(backtrack,lambdaBackTrack,Gap,Body1,Body2,deltaf,jj,imax);

                    end
                    if conv % to be sure we have converged
                        break;
                    end

                end
                if ~conv % remember it , for paper checking
                    failed = true;
                    break;
                end
                if  approachSubtype == "Augumented Lagrange" || approachSubtype == "penalty-Nitsche" || ...
                    approachSubtype == "Augumented Nitsche"

                    lambda_old = approach.lambda.meaning;
                    [~,lambda_next] = approach.AimFunction(Body1,Body2,approach);
                    approach.lambda.converge = norm(lambda_next-lambda_old)/ max(1,norm(lambda_next)) < tol;
                    approach.lambda.meaning = lambda_next;

                else
                    approach.lambda.converge = true;
                end
        end
        Body1 = SaveResults(Body1,ii,"last"); % options: "all", "last", each by (number)
        Body2 = SaveResults(Body2,ii,"last");

        if ~approach.lambda.converge % outer loop did not converge
            failed = true;
            break;
        end
        if failed
            break
        end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MeanPenetration, MaxPenetration] = Gapfunc(Body1,Body2);   % [m]

% contact force = vertical reaction at the prescribed top edge of the indenter, per unit thickness
Body1 = Elastic(Body1);                                                    % internal forces of the final state
ContactForcePerThickness = abs(sum(Body1.Fint.vec(AppliedDispVertical)))/Body1.Lz;   % [N/m]

% Hertz line contact for this force, plane stress: 1/E* = 1/E1 + 1/E2
EffectiveModulus = 1/(1/Body1.E + 1/Body2.E);
HertzHalfWidth   = sqrt(4*ContactForcePerThickness*Radius/(pi*EffectiveModulus));   % a  [m]
HertzMaxPressure = 2*ContactForcePerThickness/(pi*HertzHalfWidth);                  % p0 [Pa]
ContactCentreX   = Body1.shift.x + Body1.Lx/2;                                      % lowest point of the arc

% FE contact pressure: -sigma_yy on the top edge of the block (eta = +1) at xi = +-1/sqrt(3)
TopNodes = Body2.contact.nodalid;
PressurePositions = [];   % distance from the contact centre [m]
FEPressure = [];
for iElem = 1:size(Body2.nloc,1)
    if sum(ismember(Body2.nloc(iElem,:), TopNodes)) == 2          % element with an edge on the contact surface
        [X,U] = GetCoorDisp(iElem,Body2.nloc,Body2.P0,Body2.u);
        for xiPoint = [-1 1]/sqrt(3)
            Sigma = Sigma_2412(Body2.E,Body2.nu,U,X,xiPoint,1,Body2.StressFunction);   % eta = +1 is the top edge
            CurrentPosition = Nm_2412(xiPoint,1)*(X+U);
            PressurePositions(end+1) = CurrentPosition(1) - ContactCentreX;    %#ok<SAGROW>
            FEPressure(end+1) = -Sigma(2,2);                                     %#ok<SAGROW>
        end
    end
end
HertzPressure = HertzMaxPressure*sqrt(max(0, 1 - (PressurePositions/HertzHalfWidth).^2));

MaxPressureError  = max(abs(FEPressure - HertzPressure))/HertzMaxPressure;   % profile error / p0
PeakPressureError = max(FEPressure)/HertzMaxPressure - 1;                     % signed


% is the Hertz comparison valid for this run?
TopEdgeX = Body2.P0(Body2.contact.nodalid,1) - ContactCentreX;            % block top nodes, measured from the centre
ElementsOverHalfWidth = nnz(TopEdgeX > 0 & TopEdgeX <= HertzHalfWidth);   % block elements between the centre and x = a
fprintf('\nHertz: P = %.4g N/m, a = %.4g mm, p0 = %.4g MPa, a/R = %.3f, block elements over a = %d\n', ...
        ContactForcePerThickness, HertzHalfWidth*1e3, HertzMaxPressure*1e-6, ...
        HertzHalfWidth/Radius, ElementsOverHalfWidth);
if HertzHalfWidth/Radius > 0.2
    warning('Hertz: a/R = %.2f > 0.2, the contact is not small compared with the radius', HertzHalfWidth/Radius);
end
if HertzHalfWidth > Body1.Lx/2
    warning('Hertz: a = %.2f mm exceeds the half-width of the indenter (%.2f mm)', HertzHalfWidth*1e3, Body1.Lx/2*1e3);
end
if ElementsOverHalfWidth < 4
    warning('Hertz: only %d block elements over a, the pressure profile is not resolved', ElementsOverHalfWidth);
end

fprintf('max |p_FE - p_Hertz| / p0 = %.3g, peak pressure error = %.3g, max penetration / h = %.3g\n', ...
        MaxPressureError, PeakPressureError, MaxPenetration/approach.h);

AccuracyError = Inf;                                        % failed run = never counted as accurate
if ~failed
    AccuracyError = MaxPenetration/approach.h;
end
AccuracyType = "max penetration / h";

AlphaEffective = alpha;                                     % penalty-Nitsche raises its own alpha
if approach.Name == "penalty-Nitsche"
    AlphaEffective = max(alpha, max(alpha*approach.lambda.meaning/approach.penalty));
end

% same first 10 fields as TaylorPatchTest, Hertz-specific ones at the end (the sweep's catch struct has the same)
Result = struct('NewtonIterations',     Niter, ...
                'AccuracyError',        AccuracyError, ...
                'AccuracyType',         AccuracyType, ...
                'MaxPressureError',     MaxPressureError, ...
                'MeanPenetration_h',    MeanPenetration/approach.h, ...
                'MaxPenetration_h',     MaxPenetration/approach.h, ...
                'MaxPenetration_mm',    MaxPenetration*1e3, ...
                'Failed',               failed, ...
                'CPUTime_s',            titertot, ...
                'AlphaEffective',       AlphaEffective, ...
                'PeakPressureError',    PeakPressureError, ...
                'ContactForce_N_per_m', ContactForcePerThickness, ...
                'HalfWidth_to_Radius',  HertzHalfWidth/Radius);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %##################### Post-Processing ######################
if ~SWEEP
    PostProcess(Body1, Body2, ShowVisualization, WhatoToShow, ShowNodeNumbers, approach, ContactPointfunc, Gapfunc,Solution, PointsofInterest);

    [SortedPositions, Order] = sort(PressurePositions);
    HertzX = linspace(-1, 1, 201);
    figure('WindowStyle', 'normal');
    plot(SortedPositions/HertzHalfWidth, FEPressure(Order)/HertzMaxPressure, 'o-'); hold on
    plot(HertzX, sqrt(1 - HertzX.^2), 'k-', 'LineWidth', 1.2);
    xlim([-2 2]); grid on
    xlabel('x / a'); ylabel('p / p_0');
    legend(approach.Name, 'Hertz');

    
    

    % contact pressure of the method at the contact points (as plotted by Hu / Mlika) and Hertz
    [ContactPointPositions, ContactTraction] = ContactPressureAtPoints(Body1, Body2, approach);
    ContactPointX = ContactPointPositions(:,1) - ContactCentreX;
    HertzTraction = HertzMaxPressure*sqrt(max(0, 1 - (ContactPointX/HertzHalfWidth).^2));
    fprintf('contact pressure: max |t_N - p_Hertz| / p0 = %.3g, peak error = %.3g\n', ...
            max(abs(ContactTraction - HertzTraction))/HertzMaxPressure, max(ContactTraction)/HertzMaxPressure - 1);

    figure('WindowStyle', 'normal');
    stem(ContactPointX*1e3, ContactTraction*1e-9, 'r', 'filled', 'MarkerSize', 3); hold on
    HertzX = linspace(-HertzHalfWidth, HertzHalfWidth, 201);
    plot(HertzX*1e3, HertzMaxPressure*sqrt(1 - (HertzX/HertzHalfWidth).^2)*1e-9, 'b-', 'LineWidth', 1.2);
    xlim([-1.6 1.6]*HertzHalfWidth*1e3); grid on
    xlabel('x (mm)'); ylabel('Contact pressure (GPa)');
    legend(approach.Name + ", \alpha = " + alpha, 'Hertz');

    % keep the profile for the paper figure
    save("HertzProfile_" + strrep(approach.Name," ","_") + "_alpha" + alpha + "_theta" + approach.theta + ".mat", ...
         'ContactPointX', 'ContactTraction', 'PressurePositions', 'FEPressure', 'HertzHalfWidth', 'HertzMaxPressure');
end