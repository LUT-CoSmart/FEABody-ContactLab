clc, clear, close all
% Run from SweepAlpha.m: SWEEP = true, method and alpha come from the sweep.
if ~exist('SWEEP','var') || ~SWEEP
    clc, clear, close all
    SWEEP = false;
end
clear Body1 Body2 Body2Left Body2Right approach backtrack   % never reuse bodies from a previous run
format long;
addpath(genpath(pwd));

p = 10e3;   
%########## Bodies Data : Element positioning (from (0.0) coord. ) ########
% Body 1
Body1.Lx = 0.05;
Body1.Ly = 0.01;
Body1.Lz = 0.1;
Body1.E=2e11;
Body1.nu=0.3;
Body1 = SetMaterial(Body1, "KS"); 
Body1.shift.y = Body1.Ly - 1e-10;   % tiny initial overlap -> contact is active at iteration 1

% Body 2
Body2.Lx = 0.05;
Body2.Ly = 0.01;
Body2.Lz = 0.1;
Body2.E=2e11;
Body2.nu=0.3;
Body2 = SetMaterial(Body2, "KS"); 

%#################### Mesh #########################################
dx1 = 16;
dy1 = 4;

dx2 = 10;
dy2 = 8;

Body1.nElems.x = dx1;
Body1.nElems.y = dy1;

Body2.nElems.x = dx2;
Body2.nElems.y = dy2;

Body1 = CreateFEMesh(Body1);
Body2 = CreateFEMesh(Body2);

%#################### BC ###########################################
Body1.loc.x = 0;
Body1.loc.y = 'all';
Body1 = CreateBC(Body1,'ux');


Body2.loc.x = 'all';
Body2.loc.y = 0;
Body2 = CreateBC(Body2,'uy');

Body2.loc.x = 0;
Body2.loc.y = 'all';
Body2 = CreateBC(Body2,'ux');

%##################### Loadings ####################################
Body1.Fext.x = 0;
Body1.Fext.y = -p*Body1.Lx*Body1.Lz;
Body1.Fext.loc.x = 'all';
Body1.Fext.loc.y = Body1.Ly;


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
    approachSubtype = "Nitsche";
    alpha = 0.10;
    %##################### Newton iter. parameters ######################
    imax = 50; 
    tol=1e-3;   
    LoadType = "cubic"; % Update forces, supported loading types: linear, exponential, quadratic, cubic;
    steps= 20;
    Solution = "Newton-Rapson"; % Options: Newton-Rapson, Newton-Broyden
    PenetrationTol = 1e-7;
    ShowVisualization = true;
    ShowNodeNumbers = true;
    WhatoToShow = "sigma_yy"; % options: "ux", "uy", "u_total", "sigma_xx", "sigma_yy", "sigma_xy", "sigma_VM" 
end

PointsofInterest.Name = "Gauss";
PointsofInterest.n = 2;

[ContactPointfunc, Gapfunc, GapfuncPairs]  = ContactPointSetting(PointsofInterest,approachBasis);
Perturbation = "automatic"; % Options: "automatic", "incremental"

[approach,backtrack] = ApproachSettings(approachBasis,approachSubtype,ContactPointfunc,...
                       GapfuncPairs,Perturbation,PointsofInterest,Body1,Body2,alpha,PenetrationTol);
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

        Body1 = CreateFext(ii,steps,Body1,LoadType);
        Body2 = CreateFext(ii,steps,Body2,LoadType);


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

% Exact small-strain, plane-stress solution of this patch test (same E and nu in both blocks):
%   u_x = nu*p/E*x,   u_y = -p/E*y      (u_x = 0 along the left edge x = 0 of both blocks)
StrainX = Body1.nu*p/Body1.E;
StrainY = -p/Body1.E;
ExactU1 = [StrainX*Body1.P0(:,1), StrainY*Body1.P0(:,2)]';   % 2 x number of nodes
ExactU2 = [StrainX*Body2.P0(:,1), StrainY*Body2.P0(:,2)]';
NodalError   = [vecnorm(reshape(Body1.u,2,[]) - ExactU1), vecnorm(reshape(Body2.u,2,[]) - ExactU2)];
ErrorValue   = max(NodalError) / max([vecnorm(ExactU1), vecnorm(ExactU2)]);
AccuracyType = "max nodal displacement error / max exact displacement";

% Max contact pressure relative error on the lower contact surface (as in Mlika et al. 2017)
TopNodes = Body2.contact.nodalid;
MaxPressureError = 0;
for iElem = 1:size(Body2.nloc,1)
    if sum(ismember(Body2.nloc(iElem,:), TopNodes)) == 2          % element with an edge on the contact surface
        [X,U] = GetCoorDisp(iElem,Body2.nloc,Body2.P0,Body2.u);
        for xiPoint = [-1 1]/sqrt(3)
            Sigma = Sigma_2412(Body2.E,Body2.nu,U,X,xiPoint,1,Body2.StressFunction);   % eta = +1 is the top edge
            MaxPressureError = max(MaxPressureError, abs(-Sigma(2,2)/p - 1));
        end
    end
end
MaxPressureError
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
                'MaxPressureError',  MaxPressureError, ...
                'MeanPenetration_h', MeanPenetration/approach.h, ...
                'MaxPenetration_h',  MaxPenetration/approach.h, ...
                'MaxPenetration_mm', MaxPenetration*1e3, ...
                'Failed',            failed, ...
                'CPUTime_s',         titertot, ...
                'AlphaEffective',    AlphaEffective);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %##################### Post-Processing ######################
if ~SWEEP
    PostProcess(Body1, Body2, ShowVisualization, WhatoToShow, ShowNodeNumbers, approach, ContactPointfunc, Gapfunc,Solution, PointsofInterest);
end