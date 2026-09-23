clc, clear, close all
format long;
addpath(genpath(pwd));

%########## Element positioning (from (0.0) coord. ) ###########
% Body 1
Body1.Lx = 2;
Body1.Ly = 0.5;
Body1.Lz = 0.1;
Body1.E=2.07e11;
Body1.nu=0.3;
Body1.shift.x = 0;
Body1.shift.y = 0;

% Body 2
Body2.Lx = 2;
Body2.Ly = 0.5;
Body2.Lz = 0.1;
Body2.E=2.07e11;
Body2.nu=0.3;
Body2.shift.x = 0;
Body2.shift.y = -Body2.Ly;

%#################### Mesh #########################################
dx1 = 8;
dy1 = 3;

dx2 = 10;
dy2 = 2;

Body1.nElems.x = dx1;
Body1.nElems.y = dy1;

Body2.nElems.x = dx2;
Body2.nElems.y = dy2;

Body1 = CreateFEMesh(Body1);
Body2 = CreateFEMesh(Body2);

%#################### BC ###########################################
Body1.loc.x = 0; 
Body1.loc.y = 'all';  % Location of nodes along the axis or 'all' can be an option

Body2.loc.x = 0; 
Body2.loc.y = 'all'; 

BCtype = 'all'; % options: 'all', 'ux', 'uy' (P.S. possible to combine calling second time CreateBC)
Body1 = CreateBC(Body1,BCtype);
Body2 = CreateBC(Body2,BCtype); 

%##################### Loadings ######################
% local positions (assuming each body inside the square [(0.0)-(Lx,Ly)]
Body1.Fext.x = 0; 
Body1.Fext.y = -62.5*10^(6)*2;

Body1.Fext.loc.x = Body1.Lx;
Body1.Fext.loc.y = 'all';

Body2.Fext.y = 0; 
Body2.Fext.x = 0;

Body2.Fext.loc.x = Body2.Lx;
Body2.Fext.loc.y = 'all';

%##################### Identification of possble contact surfaces #########
% local positions (assuming all bodies in (0,0) )
Body1.contact.loc.x = 'all';
Body1.contact.loc.y = 0;   
Body1.contact.nodalid = FindGlobNodalID(Body1.P0,Body1.contact.loc,Body1.shift);

Body2.contact.loc.x = 'all';
Body2.contact.loc.y = Body2.Ly;  
Body2.contact.nodalid = FindGlobNodalID(Body2.P0,Body2.contact.loc,Body2.shift);

%##################### Contact ############################
approachBasis = "Penalty";  % Options: None, Penalty, Lagrange
approachSubtype = "Penalty";   % Subtypes
                                         % Penalty: Penalty, Augumented Lagrange,
                                         %          Nitsche, penalty-Nitsche,  Augmented Nitsche
                                         % Lagrange: Lagrange, perturbed Lagrange

PointsofInterest.Name = "nodes"; % options: "nodes", "Gauss", "LinSpace" 
PointsofInterest.n = 1; % number of points per segment (Gauss & LinSpace points)

% ==============================================================================================================
% For Lagrange-based methods, a large number of contact points can lead to an overconstrained solution.
% The reasoning is as follows. The points can be projected over several elements of the target body, which shape functions are linear.
% That makes several consequenced elements be in one line. The same reason is to have coarser mesh for the contact body than the target one.  
% ==============================================================================================================

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
type = "cubic"; % Update forces, supported loading types: linear, exponential, quadratic, cubic;
steps= 5;
Solution = "Newton-Rapson"; % Options: Newton-Rapson, Newton-Broyden
%#################### Processing ######################

total_steps = 0;
titertot=0;  

for ii = 1:steps

        approach.lambda.step = 0;
        approach.lambda.converge = false;  
        
        Body1 = CreateFext(ii,steps,Body1,type);
        Body2 = CreateFext(ii,steps,Body2,type);
        
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
WhatoToShow = "u_total"; % options: "ux", "uy", "u_total", "sigma_xx", "sigma_yy", "sigma_xy" 
PostProcess(Body1, Body2, ShowVisualization, WhatoToShow, ShowNodeNumbers, approach, ContactPointfunc, Gapfunc,Solution, PointsofInterest);

