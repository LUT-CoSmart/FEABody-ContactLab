clc, clear, close all
if ~exist('SWEEP','var')|| ~SWEEP % we have apper file to deal with    
    clc, clear, close all    
    SWEEP = false;
else
    clear Body1 Body2 approach backtrack Result Assemblance BestStep Contact
end

format long;
addpath(genpath(pwd));

p = 2e9;
W = 6;

%########## Bodies Data : Element positioning (from (0.0) coord. ) ########
% Body 1
Body1.Lx = 1.2;
Body1.Ly = 4;
Body1.Lz = 0.1;
Body1.E=68.96*10^(8+6);
Body1.nu=0.32;
Radius = 0.75;
Body1 = SetMaterial(Body1, "Neo"); 

% Body 2
Body2.Lx = 12;
Body2.Ly = 5;
Body2.Lz = 0.1;
Body2.E=68.96*10^(7+6);
Body2.nu=0.32;
Body2 = SetMaterial(Body2, "Neo"); 


Body1.shift.x = W;
Body1.shift.y = Body2.Ly - 1e-11;   % tiny initial overlap -> contact is active at iteration 1

%#################### Mesh #########################################
dx1 = 7;
dy1 = 2;

dx2 = 8;
dy2 = 8;

Body1.nElems.x = dx1;
Body1.nElems.y = dy1;

Body2.nElems.x = dx2;
Body2.nElems.y = dy2;

Body1 = CreateFEMesh(Body1);
Body2 = CreateFEMesh(Body2);


% Visualization(Body1);
% Visualization(Body2);

 
% %#################### BC ###########################################
Body1.loc.x = 'all';
Body1.loc.y = Body1.Ly;
Body1 = CreateBC(Body1,'all');

Body2.loc.x = 'all';
Body2.loc.y = 0;
Body2 = CreateBC(Body2,'all');

%##################### Prescribed displacement of the die top (instead of a force) ##############
VerticalDisplacement = 0.01;
AppliedDisp.x = 'all';
AppliedDisp.y = Body1.Ly;
AppliedDispNodes = FindGlobNodalID(Body1.P0, AppliedDisp, Body1.shift);
AppliedDispVerticalDofs = xlocChosen(Body1.DofsAtNode, AppliedDispNodes, 2);

%##################### Contact surfaces #############################
Body1.contact.loc.x = 'all';
Body1.contact.loc.y = 0;
Body1.contact.nodalid = FindGlobNodalID(Body1.P0,Body1.contact.loc,Body1.shift);

Body2.contact.loc.x = 'all';
Body2.contact.loc.y = Body2.Ly;
Body2.contact.nodalid = FindGlobNodalID(Body2.P0,Body2.contact.loc,Body2.shift);

%##################### Making bottom surface round #############################
Body1 = RoundingBottom(Body1, Radius);

%##################### Contact ######################################

% Subtypes
% Penalty: Penalty, Augumented Lagrange,
%          Nitsche, penalty-Nitsche,  Augumented Nitsche
% Lagrange: Lagrange, perturbed Lagrange
if ~SWEEP 
    approachBasis = "Penalty";
    approachSubtype = "Penalty";
    alpha = 0.10;
    %##################### Newton iter. parameters ######################
    imax = 50; 
    tol=1e-3;   
    LoadType = "linear"; % Update forces, supported loading types: linear, exponential, quadratic, cubic;
    steps= 10;
    Solution = "Newton-Rapson"; % Options: Newton-Rapson, Newton-Broyden
    PenetrationTol = 1e-7;
    ShowVisualization = true;
    ShowNodeNumbers = true;
    WhatoToShow = "sigma_yy"; % options: "ux", "uy", "u_total", "sigma_xx", "sigma_yy", "sigma_xy", "sigma_VM" 
end

PointsofInterest.Name = "Nodes";
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

        Body1 = CreateFext(ii,steps,Body1,LoadType);  % these give zeros, but necessary
        Body2 = CreateFext(ii,steps,Body2,LoadType);  % these give zeros, but necessary
        Body1.u(AppliedDispVerticalDofs) = -VerticalDisplacement*ii/steps;   % constrained DOFs: Newton never changes them, so this value is kept

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
% [MeanPenetration, MaxPenetration] = Gapfunc(Body1,Body2);   % [m]
% 
% 
% AlphaEffective = alpha;                                     % penalty-Nitsche raises its own alpha
% if approach.Name == "penalty-Nitsche"
%     AlphaEffective = max(alpha, max(alpha*approach.lambda.meaning/approach.penalty));
% end
% 
% Result = struct('NewtonIterations',  Niter, ...
%                 'AccuracyError',     AccuracyError, ...
%                 'AccuracyType',      AccuracyType, ...
%                 'MeanPenetration_h', MeanPenetration/approach.h, ...
%                 'MaxPenetration_h',  MaxPenetration/approach.h, ...
%                 'Failed',            failed, ...
%                 'CPUTime_s',         titertot, ...
%                 'AlphaEffective',    AlphaEffective);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %##################### Post-Processing ######################
if ~SWEEP
    PostProcess(Body1, Body2, ShowVisualization, WhatoToShow, ShowNodeNumbers, approach, ContactPointfunc, Gapfunc,Solution, PointsofInterest);
end