%% From Kim J. H. et Al., 1976
clc, clear, close all
if ~exist('SWEEP','var')|| ~SWEEP % we have apper file to deal with    
    clc, clear, close all    
    SWEEP = false;
else
    clear Body1 Body2 approach backtrack Result Assemblance BestStep Contact
end

format long;
addpath(genpath(pwd));

p = 2e7;    % pressure [Pa]. Small (p/E = 1e-4) so the exact small-strain solution below is valid
W = 5;

%########## Bodies Data : Element positioning (from (0.0) coord. ) ########
% Body 1
Body1.Lx = 10;
Body1.Ly = 5;
Body1.Lz = 0.1;
Body1.E=2e11;
Body1.nu=0.28;
Body1 = SetMaterial(Body1, "KS"); 

% Body 2
Body2.Lx = 20;
Body2.Ly = 5;
Body2.Lz = 0.1;
Body2.E=2e11;
Body2.nu=0.28;
Body2 = SetMaterial(Body2, "KS"); 

Body1.shift.x = W;
Body1.shift.y = Body2.Ly - 1e-7;   % tiny initial overlap -> contact is active at iteration 1
%#################### Mesh #########################################
dx1 = 8;
dy1 = 4;

dx2 = 20;
dy2 = 8;

Body1.nElems.x = dx1;
Body1.nElems.y = dy1;

Body2.nElems.x = dx2;
Body2.nElems.y = dy2;

Body1 = CreateFEMesh(Body1);
Body2 = CreateFEMesh(Body2);

%#################### BC ###########################################
Body1.loc.x = Body1.Lx/2;
Body1.loc.y = 0;
Body1 = CreateBC(Body1,'ux');

% Lower rollers: vertical displacement fixed along the bottom.
Body2.loc.x = 'all';
Body2.loc.y = 0;
Body2 = CreateBC(Body2,'uy');

% Lower midpoint pin: additionally fix its horizontal displacement.
Body2.loc.x = Body2.Lx/2;
Body2.loc.y = 0;
Body2 = CreateBC(Body2,'ux');

%##################### Loadings ####################################
% Whole top edge of the upper block.
Body1.Fext.x = 0;
Body1.Fext.y = -p*Body1.Lx*Body1.Lz;
Body1.Fext.loc.x = 'all';
Body1.Fext.loc.y = Body1.Ly;

% Same pressure on the exposed top parts of the lower block (left and right of the upper block).
% Needed for the exact patch-test solution (uniform sigma_yy = -p in both bodies).
% Mesh requirement: x = W and x = W + Body1.Lx must be nodes of Body2.
LoadExposedParts = true;
if LoadExposedParts
    Body2Left = Body2;                                   % copy, used only to build the load on [0, W]
    Body2Left.Fext.x = 0;
    Body2Left.Fext.y = -p*W*Body2.Lz;
    Body2Left.Fext.loc.x = [0 W];
    Body2Left.Fext.loc.y = Body2.Ly;

    Body2Right = Body2;                                  % copy, used only to build the load on [W+Lx1, Lx2]
    Body2Right.Fext.x = 0;
    Body2Right.Fext.y = -p*(Body2.Lx - W - Body1.Lx)*Body2.Lz;
    Body2Right.Fext.loc.x = [W+Body1.Lx, Body2.Lx];
    Body2Right.Fext.loc.y = Body2.Ly;
end

%##################### Contact surfaces #############################
Body1.contact.loc.x = 'all';
Body1.contact.loc.y = 0;
Body1.contact.nodalid = FindGlobNodalID(Body1.P0,Body1.contact.loc,Body1.shift);

Body2.contact.loc.x = 'all';
Body2.contact.loc.y = Body2.Ly;
Body2.contact.nodalid = FindGlobNodalID(Body2.P0,Body2.contact.loc,Body2.shift);

%##################### Contact ######################################

% Subtypes
% Penalty: Penalty, Augumented Lagrange,
%          Nitsche, penalty-Nitsche,  Augumented Nitsche
% Lagrange: Lagrange, perturbed Lagrange
if ~SWEEP 
    approachBasis = "Penalty";
    approachSubtype = "penalty-Nitsche";
    alpha = 0.10;
    %##################### Newton iter. parameters ######################
    imax = 50; 
    tol=1e-3;   
    LoadType = "cubic"; % Update forces, supported loading types: linear, exponential, quadratic, cubic;
    steps= 1;
    Solution = "Newton-Rapson"; % Options: Newton-Rapson, Newton-Broyden
    PenetrationTol = 1e-7;
    ShowVisualization = true;
    ShowNodeNumbers = true;
    WhatoToShow = "sigma_yy"; % options: "ux", "uy", "u_total", "sigma_xx", "sigma_yy", "sigma_xy", "sigma_VM" 
end

PointsofInterest.Name = "Gauss";
PointsofInterest.n = 1;

[ContactPointfunc, Gapfunc, GapfuncPairs]  = ContactPointSetting(PointsofInterest,approachBasis);
Perturbation = "automatic"; % Options: "automatic", "incremental"

[approach,backtrack] = ApproachSettings(approachBasis,approachSubtype,ContactPointfunc,...
                       GapfuncPairs,Perturbation,PointsofInterest,Body1,Body2,alpha,PenetrationTol);
if ~SWEEP
    approach.theta = -1;
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

        Body1 = CreateFext(ii,steps,Body1,LoadType);
        if LoadExposedParts
            Body2Left  = CreateFext(ii,steps,Body2Left,LoadType);
            Body2Right = CreateFext(ii,steps,Body2Right,LoadType);
            Body2.Fext = Body2Left.Fext;                                   % keeps nodalid for SaveResults
            Body2.Fext.vec = Body2Left.Fext.vec + Body2Right.Fext.vec;     % load on both exposed parts
        else
            Body2 = CreateFext(ii,steps,Body2,LoadType);
        end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Metrics for the alpha sweep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MeanPenetration, MaxPenetration] = Gapfunc(Body1,Body2);   % [m]

if LoadExposedParts
    % Exact small-strain, plane-stress solution (both bodies have the same E and nu):
    %   u_x = nu*p/E*(x - x_pin),   u_y = -p/E*y
    % x_pin = x of the horizontally fixed node of each body (frictionless contact -> each body has its own)
    StrainX = Body1.nu*p/Body1.E;
    StrainY = -p/Body1.E;
    PinX1   = Body1.shift.x + Body1.Lx/2;
    PinX2   = Body2.shift.x + Body2.Lx/2;
    ExactU1 = [StrainX*(Body1.P0(:,1) - PinX1), StrainY*Body1.P0(:,2)]';   % 2 x number of nodes
    ExactU2 = [StrainX*(Body2.P0(:,1) - PinX2), StrainY*Body2.P0(:,2)]';
    NodalError = [vecnorm(reshape(Body1.u,2,[]) - ExactU1), vecnorm(reshape(Body2.u,2,[]) - ExactU2)];
    ErrorValue   = max(NodalError) / max([vecnorm(ExactU1), vecnorm(ExactU2)]);
    AccuracyType = "max nodal displacement error / max exact displacement";
else
    % no exact solution without the load on the exposed parts -> constraint satisfaction only
    ErrorValue   = MaxPenetration/approach.h;
    AccuracyType = "max penetration / h";
end

AccuracyError = Inf;                                        % failed run = never counted as accurate
if ~failed
    AccuracyError = ErrorValue;
end

AlphaEffective = alpha;                                     % penalty-Nitsche raises its own alpha
if approach.Name == "penalty-Nitsche"
    AlphaEffective = max(alpha, max(alpha*approach.lambda.meaning/approach.penalty));
end

Result = struct('NewtonIterations',  Niter, ...
                'AccuracyError',     AccuracyError, ...
                'AccuracyType',      AccuracyType, ...
                'MeanPenetration_h', MeanPenetration/approach.h, ...
                'MaxPenetration_h',  MaxPenetration/approach.h, ...
                'Failed',            failed, ...
                'CPUTime_s',         titertot, ...
                'AlphaEffective',    AlphaEffective);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %##################### Post-Processing ######################
if ~SWEEP
    PostProcess(Body1, Body2, ShowVisualization, WhatoToShow, ShowNodeNumbers, approach, ContactPointfunc, Gapfunc,Solution, PointsofInterest);
end