% ShowSweepResults.m - reads the results saved by SweepAlpha.m and shows the alpha regions.
% Nothing is computed again; only the saved Table is analysed.
clc, clear, close all

CaseName        = 'LayredBeams';   % 'PatchTest' or 'LayredBeams'
AccuracyTol     = 1e-3;            % accurate if max penetration / h <= AccuracyTol (same tol as in the sweep)
TipDeviationTol = 1e-2;            % used only if the sweep recorded TipDeflection

load(['Sweep_' CaseName '.mat'], 'Table');

Methods = unique(Table.Method, 'stable');
Alphas  = unique(Table.Alpha);

%% ---------------- Converged / correct per run ----------------
Table.Converged = ~Table.Failed;
Table.Correct   = Table.Converged & Table.AccuracyError <= AccuracyTol;

HasTipDeflection = ismember('TipDeflection', Table.Properties.VariableNames) && any(~isnan(Table.TipDeflection));
if HasTipDeflection
    ReferenceTipDeflection = median(Table.TipDeflection(Table.Correct));
    Table.TipDeviation     = abs(Table.TipDeflection./ReferenceTipDeflection - 1);
    Table.Correct          = Table.Correct & Table.TipDeviation <= TipDeviationTol;
end
disp(Table)

%% ---------------- Regions per method (widest correct region first) ----------------
Summary = table();
for iM = 1:numel(Methods)
    T = sortrows(Table(Table.Method == Methods(iM),:), 'Alpha');
    [ConvergedRanges, ConvergedDecades] = AlphaRanges(T.Alpha, T.Converged);
    [CorrectRanges,   CorrectDecades]   = AlphaRanges(T.Alpha, T.Correct);
    Summary = [Summary; table(Methods(iM), ConvergedRanges, ConvergedDecades, CorrectRanges, CorrectDecades, ...
               'VariableNames', {'Method','ConvergedAlphaRanges','ConvergedDecades','CorrectAlphaRanges','CorrectDecades'})];
end
Summary = sortrows(Summary, 'CorrectDecades', 'descend');
disp(Summary)
fprintf('Widest correct region:    %s (%.2f decades)\n', Summary.Method(1),   Summary.CorrectDecades(1));
fprintf('Narrowest correct region: %s (%.2f decades)\n', Summary.Method(end), Summary.CorrectDecades(end));

%% ---------------- Fair comparison: alphas where every method is correct ----------------
IsFairAlpha = arrayfun(@(a) sum(Table.Alpha == a) == numel(Methods) && all(Table.Correct(Table.Alpha == a)), Alphas);
FairAlphas  = Alphas(IsFairAlpha);
if isempty(FairAlphas)
    fprintf('No alpha where all methods give the correct solution.\n');
else
    fprintf('Alphas where all methods give the correct solution: %s\n', mat2str(FairAlphas.', 3));
end
FairColumns = intersect({'Method','Alpha','NewtonIterations','CPUTime_s','MaxPenetration_h','TipDeviation','AlphaEffective'}, ...
                        Table.Properties.VariableNames, 'stable');
FairTable   = sortrows(Table(ismember(Table.Alpha, FairAlphas), FairColumns), {'Alpha','Method'});
disp(FairTable)

function [RangesText, LongestBlockDecades] = AlphaRanges(Alpha, IsOk)
    % contiguous alpha blocks where IsOk is true, e.g. "[0.00316 0.0316] [0.316 1e+04]",
    % and the width of the longest block in decades
    Change     = diff([0; IsOk(:); 0]);
    BlockStart = find(Change == 1);
    BlockEnd   = find(Change == -1) - 1;
    RangesText = "none";
    LongestBlockDecades = 0;
    if ~isempty(BlockStart)
        RangesText = strjoin(compose("[%g %g]", Alpha(BlockStart), Alpha(BlockEnd)).', " ");
        LongestBlockDecades = max(log10(Alpha(BlockEnd)./Alpha(BlockStart)));
    end
end