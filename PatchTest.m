clc, clear, close all
format long;
addpath(genpath(pwd));

p = 1;    
W = 1;

%########## Bodies Data : Element positioning (from (0.0) coord. ) ########
% Body 1
Body1.Lx = 10;
Body1.Ly = 5;
Body1.Lz = 0.1;
Body1.E=2e11;
Body1.nu=0.28;
Body1.shift.x = W;
Body1.shift.y = 5;

% Body 2
Body2.Lx = 20;
Body2.Ly = 5;
Body2.Lz = 0.1;
Body2.E=2e11;
Body2.nu=0.28;
Body2.shift.x = 0;
Body2.shift.y = 0;

%#################### Mesh #########################################
dx1 = 4;
dy1 = 1;

dx2 = 4;
dy2 = 1;

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

% % Lower midpoint pin: additionally fix its horizontal displacement.
Body2.loc.x = Body2.Lx/2;
Body2.loc.y = 0;
Body2 = CreateBC(Body2,'ux');

%##################### Loadings ####################################
% Whole top edge of the upper block.
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
approachBasis = "Penalty";
approachSubtype = "Penalty"; 
PointsofInterest.Name = "Gauss";
PointsofInterest.n = 2;

[ContactPointfunc, Gapfunc, GapfuncPairs]  = ContactPointSetting(PointsofInterest);
Perturbation = "automatic"; % Options: "automatic", "incremental"
[approach,backtrack] = ApproachSettings(approachBasis,approachSubtype,ContactPointfunc,...
                       GapfuncPairs,Perturbation,PointsofInterest);

if approachSubtype == "Augumented Lagrange" || approachSubtype == "penalty-Nitsche" || approachSubtype == "Augumented Nitsche"  
    [InitialContactPoints,~] = ContactPointfunc(Body1);
    approach.lambda.meaning = zeros(size(InitialContactPoints,1),1);
    approach.lambda.imax = 100;
end

%##################### Newton iter. parameters ######################
imax = 20; 
tol=1e-4;   
LoadType = "cubic"; % Update forces, supported loading types: linear, exponential, quadratic, cubic;
steps= 1;
Solution = "Newton-Rapson"; % Options: Newton-Rapson, Newton-Broyden
%#################### Processing ######################

total_steps = 0;
titertot=0;  


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
                        % interaction of two bodies
                        [DofsFunction, Stiffness] = Contact(Body1,Body2,approach,jj,Solution);
        
                        Gap = Gapfunc(Body1,Body2);
                        % inner forces of the each body
                        Body1 = Elastic(Body1);
                        Body2 = Elastic(Body2);

                        [Body1, Body2, uu_bc, deltaf, lambda_next] = Assemblance(jj,Solution,Body1, Body2, DofsFunction,Stiffness,approach,lambdaBackTrack);
                      
                        titer=toc;
                        titertot=titertot+titer;
        
                        if printStatus(approachBasis, deltaf, uu_bc, tol, ii, jj, imax, steps, titertot, Gap)
                            break;  
                        end  
    
                        [Body1,Body2] = BestStep(backtrack,lambdaBackTrack,Gap,Body1,Body2,deltaf,jj,imax);
    
                    end
    
                    if jj < imax % it is needed for exit from the second loop (exit from backtrack)
                        break;
                    end

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

end
% %##################### Post-Processing ######################
ShowVisualization = true;
ShowNodeNumbers = true;
WhatoToShow = "sigma_yy"; % options: "ux", "uy", "u_total", "sigma_xx", "sigma_yy", "sigma_xy" 
PostProcess(Body1, Body2, ShowVisualization, WhatoToShow, ShowNodeNumbers, approach, ContactPointfunc, Gapfunc,Solution, PointsofInterest);


axis equal

% Full-load check: -4e9 N for p = 2e9 Pa and thickness = 0.1 m.
Fy1 = sum(Body1.Fext.vec(2:2:end));
Fy2 = sum(Body2.Fext.vec(2:2:end));
fprintf('Applied Fy: upper = %.6e N, lower = %.6e N, total = %.6e N\n', Fy1,Fy2,Fy1 + Fy2);

