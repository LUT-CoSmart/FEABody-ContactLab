clc, clear, close all
addpath(genpath(pwd));
SWEEP = true; % used in other prgrams, when called

if SWEEP
    approachBasis = "Penalty";  % Options: None, Penalty, Lagrange 
    imax = 50;
    tol = 1e-4;              % Newton and lambda-loop tolerance (all methods)
    PenetrationTol = 1e-5;   % accuracy target: max penetration / h
    LoadType = "linear";
    steps = 10;
    Solution = "Newton-Rapson";
    BackTrack = false;
    theta = -1;
end

CaseName    = 'HertzProblem'; % 'PatchTest' or 'LayredBeams', 'TaylorPatchTest', 'HertzProblem'
Methods     = ["penalty-Nitsche"]; % "Augumented Nitsche"
if theta == 1
     Alphas      = 10.^(0:0.5:2);
     % Alphas      = 10.^[0 0.5];
elseif theta == -1
    Alphas      = 10.^(-2:0.5:2);
end
AccuracyTol = PenetrationTol;


Table = table();
for iM = 1:numel(Methods)
    for iA = 1:numel(Alphas)
        approachSubtype = Methods(iM);   
        alpha = Alphas(iA);
        [approachSubtype alpha]
        
        try            
            run([CaseName '.m']);
        catch ME
            fprintf(2,'%s, alpha = %g crashed: %s\n', Methods(iM), alpha, ME.message);
            Result = struct('NewtonIterations',NaN,'AccuracyError',Inf,'AccuracyType',"crashed", ...
                            'MaxPressureError',NaN,'MeanPenetration_h',NaN,'MaxPenetration_h',NaN, ...
                            'MaxPenetration_mm',NaN,'Failed',true,'CPUTime_s',NaN,'AlphaEffective',alpha, ...
                            'PeakPressureError',NaN,'ContactForce_N_per_m',NaN,'HalfWidth_to_Radius',NaN);
        end
        Row   = [table(Methods(iM), alpha, 'VariableNames', {'Method','Alpha'}), struct2table(Result)];
        Table = [Table; Row];
    end
end
save(['Sweep_' CaseName '.mat'],'Table');
disp(Table)

Summary = table();
for iM = 1:numel(Methods)
    T  = Table(Table.Method == Methods(iM),:);
    ok = ~T.Failed & T.AccuracyError <= AccuracyTol;
    d  = diff([0; ok; 0]);
    BlockStart = find(d == 1);
    BlockEnd   = find(d == -1) - 1;
    if isempty(BlockStart)
        AlphaMin = NaN;  AlphaMax = NaN;  MinNewtonIterations = Inf;
    else
        [~,k] = max(BlockEnd - BlockStart);
        AlphaMin = T.Alpha(BlockStart(k));
        AlphaMax = T.Alpha(BlockEnd(k));
        MinNewtonIterations = min(T.NewtonIterations(ok));
    end
    Summary = [Summary; table(Methods(iM), AlphaMin, AlphaMax, MinNewtonIterations, ...
               'VariableNames', {'Method','AlphaMin','AlphaMax','MinNewtonIterations'})];
end
disp(Summary)

clear SWEEP   % so the next stand-alone run of LayredBeams/PatchTest is not treated as a sweep