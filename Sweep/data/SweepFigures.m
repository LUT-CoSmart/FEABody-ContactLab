% SweepFigures.m - figures and LaTeX tables for the alpha sweep of the layered beams.
% Reads the tables saved by SweepAlpha_run.m (Sweep/data/*.mat); nothing is computed again.
% Correctness criterion as in SweepAlpha_call.m: converged and max penetration / h <= tolerance.
%
% Output (folder Sweep/results):
%   sweep_penetration.pdf/.png          max (top) and mean (bottom) penetration / h vs alpha
%   sweep_newton_iterations.pdf/.png    total Newton iterations vs alpha
%   sweep_alpha_window.pdf/.png         alpha range meeting each tolerance
%   sweep_heatmap_penetration.pdf/.png  method x alpha, max penetration / h
%   sweep_heatmap_iterations.pdf/.png   method x alpha, total Newton iterations
%   sweep_tables.tex                    summary table + penetration and iteration data tables
clc, clear, close all

ThisFolder   = fileparts(mfilename('fullpath'));
DataFolder   = fullfile(ThisFolder, '..', 'data');
OutputFolder = fullfile(ThisFolder, '..', 'results');
if ~isfolder(OutputFolder)
    mkdir(OutputFolder);
end

% file names -> theta (sym = +1, unsym = -1) and PenetrationTol (only penalty-Nitsche uses it)
SweepFiles = table( ...
    ["Sweep_LayredBeams_Penet_1m3_unsym.mat"; "Sweep_LayredBeams_Penet_1em5_unsym.mat"; ...
     "Sweep_LayredBeams_Penet_1em3_sym.mat";  "Sweep_LayredBeams_Penet_1em5_sym.mat"], ...
    [-1; -1; +1; +1], ...
    [1e-3; 1e-5; 1e-3; 1e-5], ...
    'VariableNames', {'FileName', 'Theta', 'PenetrationTol'});

% methods in the paper's order: name in the sweep tables, name in figures/tables, colour, marker
Methods = table( ...
    ["Penalty"; "Nitsche"; "Augumented Nitsche"; "penalty-Nitsche"; "Augumented Lagrange"], ...
    ["Penalty"; "Nitsche"; "Augmented Nitsche"; "Penalty--Nitsche"; "Augmented Lagrange"], ...
    ["#2a78d6"; "#eb6834"; "#1baf7a"; "#eda100"; "#e87ba4"], ...
    ['o'; 's'; 'd'; '^'; 'v'], ...
    'VariableNames', {'SweepName', 'DisplayName', 'Colour', 'Marker'});

AllRuns = LoadSweepRuns(DataFolder, SweepFiles);

Figures = { ...
    "sweep_penetration",         PlotQuantityVsAlpha(AllRuns, Methods, "penetration"); ...
    "sweep_newton_iterations",   PlotQuantityVsAlpha(AllRuns, Methods, "iterations"); ...
    "sweep_alpha_window",        PlotAlphaWindow(AllRuns, Methods); ...
    "sweep_heatmap_penetration", PlotSweepHeatmap(AllRuns, Methods, "penetration"); ...
    "sweep_heatmap_iterations",  PlotSweepHeatmap(AllRuns, Methods, "iterations")};
for iFigure = 1:size(Figures, 1)
    FileBase = fullfile(OutputFolder, Figures{iFigure, 1});
    exportgraphics(Figures{iFigure, 2}, FileBase + ".pdf", 'ContentType', 'vector');
    exportgraphics(Figures{iFigure, 2}, FileBase + ".png", 'Resolution', 300);
end

WriteSweepLatexTables(AllRuns, Methods, fullfile(OutputFolder, 'sweep_tables.tex'));
fprintf('Figures and tables written to %s\n', OutputFolder);


%% ======================================================================
%  Data
%  ======================================================================
function AllRuns = LoadSweepRuns(DataFolder, SweepFiles)
    AllRuns = table();
    for iFile = 1:height(SweepFiles)
        Loaded = load(fullfile(DataFolder, SweepFiles.FileName(iFile)), 'Table');
        Runs = Loaded.Table;
        Runs.Theta          = repmat(SweepFiles.Theta(iFile), height(Runs), 1);
        Runs.PenetrationTol = repmat(SweepFiles.PenetrationTol(iFile), height(Runs), 1);
        AllRuns = [AllRuns; Runs]; %#ok<AGROW>
    end
    AllRuns.Method = string(AllRuns.Method);
    AllRuns.Failed = logical(AllRuns.Failed);
end

function Rows = SweepRows(AllRuns, Theta, SweepName, PenetrationTol)
    % runs of one method; methods other than penalty-Nitsche do not use PenetrationTol,
    % so they are taken from the 1e-3 files (the 1e-5 files hold identical runs)
    if SweepName ~= "penalty-Nitsche"
        PenetrationTol = 1e-3;
    end
    IsSelected = AllRuns.Theta == Theta & AllRuns.PenetrationTol == PenetrationTol & AllRuns.Method == SweepName;
    Rows = sortrows(AllRuns(IsSelected, :), 'Alpha');
end

function Run = RunAtAlpha(Rows, Alpha)
    Run = Rows(abs(Rows.Alpha - Alpha) <= 1e-9 * Alpha, :);
end

function Curves = CurveList(Methods)
    % one entry per plotted curve / table column; penalty-Nitsche appears once per tolerance
    Curves = struct('SweepName', {}, 'DisplayName', {}, 'Colour', {}, ...
                    'Marker', {}, 'PenetrationTol', {}, 'LineStyle', {}, 'HollowMarker', {});
    for iMethod = 1:height(Methods)
        IsPenaltyNitsche = Methods.SweepName(iMethod) == "penalty-Nitsche";
        if IsPenaltyNitsche
            Tolerances = [1e-3, 1e-5];
        else
            Tolerances = 1e-3;
        end
        for Tolerance = Tolerances
            DisplayName = Methods.DisplayName(iMethod);
            if IsPenaltyNitsche
                DisplayName = DisplayName + ", $\varepsilon=10^{" + string(round(log10(Tolerance))) + "}$";
            end
            IsSecondTolerance = IsPenaltyNitsche && Tolerance < 1e-3;
            LineStyle = '-';
            if IsSecondTolerance
                LineStyle = '--';
            end
            Curves(end + 1) = struct('SweepName', Methods.SweepName(iMethod), ... %#ok<AGROW>
                'DisplayName', DisplayName, ...
                'Colour', HexToRgb(Methods.Colour(iMethod)), 'Marker', Methods.Marker(iMethod), ...
                'PenetrationTol', Tolerance, 'LineStyle', LineStyle, 'HollowMarker', IsSecondTolerance);
        end
    end
end


%% ======================================================================
%  Figures
%  ======================================================================
function Fig = PlotQuantityVsAlpha(AllRuns, Methods, Quantity)
    % Quantity: "penetration" (rows: max and mean penetration / h) or "iterations" (total Newton iterations)
    Colours = FigureColours();
    AlphaLimits = [6e-3, 1.7e3];
    switch Quantity
        case "penetration"
            Columns         = {'MaxPenetration_h', 'MeanPenetration_h'};   % mean = area-weighted over penetrating points
            YLabels         = {'max penetration $/\,h$', 'mean penetration $/\,h$'};
            FigureHeight    = 16.5;
            YLimits         = [4e-16, 2.5];
            FloorThreshold  = 1e-14;   % values below (incl. exact 0) are drawn on the floor row
            FloorRow        = 1e-15;
            NotConvergedRow = 1;
            YTicks          = [1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1];
            YTickLabels     = {'$<10^{-14}$', '$10^{-12}$', '$10^{-9}$', '$10^{-6}$', '$10^{-3}$', 'n.c.'};
            Separators      = [3e-15, 0.25];
        case "iterations"
            Columns         = {'NewtonIterations'};
            YLabels         = {'total Newton iterations (40 load steps)'};
            FigureHeight    = 9.5;
            YLimits         = [100, 1750];
            FloorThreshold  = -Inf;
            FloorRow        = NaN;
            NotConvergedRow = 1500;
            YTicks          = [100, 150, 200, 300, 500, 1000, 1500];
            YTickLabels     = {'100', '150', '200', '300', '500', '1000', 'n.c.'};
            Separators      = 1250;
    end
    Curves = CurveList(Methods);
    Fig = NewPaperFigure(17.8, FigureHeight);
    Layout = tiledlayout(Fig, numel(Columns), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    LegendHandles = gobjects(0);
    PanelLetters = 'abcd';
    for iRow = 1:numel(Columns)
        Column = Columns{iRow};
        for Theta = [-1, +1]
            Ax = nexttile(Layout);
            hold(Ax, 'on');
            if Theta == +1
                patch(Ax, [AlphaLimits(1), 0.8, 0.8, AlphaLimits(1)], YLimits([1, 1, 2, 2]), Colours.NotRunFill, ...
                      'EdgeColor', 'none', 'HandleVisibility', 'off');
                text(Ax, 0.09, sqrt(YLimits(1) * YLimits(2)), 'not run', 'Rotation', 90, ...
                     'HorizontalAlignment', 'center', 'Color', Colours.InkSecondary, 'FontSize', 8);
            end
            for Separator = Separators
                yline(Ax, Separator, '-', 'Color', Colours.Grid, 'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
            if Quantity == "penetration"
                yline(Ax, 1e-3, ':', '$\varepsilon=10^{-3}$', 'Interpreter', 'latex', 'Color', Colours.InkSecondary, ...
                      'LineWidth', 0.7, 'FontSize', 7, 'LabelHorizontalAlignment', 'right', ...
                      'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off');
                yline(Ax, 1e-5, ':', '$\varepsilon=10^{-5}$', 'Interpreter', 'latex', 'Color', Colours.InkSecondary, ...
                      'LineWidth', 0.7, 'FontSize', 7, 'LabelHorizontalAlignment', 'right', ...
                      'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
            end
            for Curve = Curves
                Rows = SweepRows(AllRuns, Theta, Curve.SweepName, Curve.PenetrationTol);
                if isempty(Rows)
                    continue
                end
                Alpha  = Rows.Alpha;
                Values = Rows.(Column);
                Failed = Rows.Failed;
                Values(Values < FloorThreshold) = FloorRow;
                MarkerFace = Curve.Colour;
                if Curve.HollowMarker
                    MarkerFace = [1, 1, 1];
                end
                CurveHandle = plot(Ax, Alpha(~Failed), Values(~Failed), 'LineStyle', Curve.LineStyle, ...
                    'Color', Curve.Colour, 'LineWidth', 1.2, 'Marker', Curve.Marker, 'MarkerSize', 5, ...
                    'MarkerFaceColor', MarkerFace, 'MarkerEdgeColor', Curve.Colour, 'DisplayName', Curve.DisplayName);
                if any(Failed)
                    plot(Ax, Alpha(Failed) * FailedMarkerShift(Curve.SweepName), repmat(NotConvergedRow, nnz(Failed), 1), ...
                         'x', 'Color', Curve.Colour, 'MarkerSize', 7, 'LineWidth', 1.4, 'Clipping', 'off', ...
                         'HandleVisibility', 'off');
                end
                if Theta == -1 && iRow == 1
                    LegendHandles(end + 1) = CurveHandle; %#ok<AGROW>
                end
            end
            set(Ax, 'XScale', 'log', 'YScale', 'log', 'XLim', AlphaLimits, 'YLim', YLimits, ...
                'XTick', 10.^(-2:3), 'XTickLabel', cellstr(compose("$10^{%d}$", (-2:3)')), ...
                'YTick', YTicks, 'YTickLabel', YTickLabels, 'YMinorTick', 'off');
            StyleAxes(Ax, Colours);
            Ax.XGrid = 'on';
            if Quantity == "iterations"
                Ax.YGrid = 'on';
            end
            if iRow == numel(Columns)
                xlabel(Ax, 'penalty parameter $\alpha$', 'Color', Colours.Ink);
            end
            if Theta == -1
                ylabel(Ax, YLabels{iRow}, 'Color', Colours.Ink);
            else
                Ax.YTickLabel = {};
            end
            PanelLetter = PanelLetters(2 * (iRow - 1) + 1 + (Theta == +1));
            title(Ax, "(" + PanelLetter + ")" + extractAfter(PanelTitle(Theta), 3), 'FontWeight', 'normal', ...
                  'FontSize', 9, 'Color', Colours.Ink);
            Ax.TitleHorizontalAlignment = 'left';
            if Theta == -1 && iRow == 1
                LegendAxes = Ax;
            end
        end
    end
    LegendHandles(end + 1) = plot(LegendAxes, NaN, NaN, 'x', 'Color', Colours.InkSecondary, 'MarkerSize', 7, ...
                                  'LineWidth', 1.4, 'DisplayName', 'not converged');
    Lgd = legend(LegendAxes, LegendHandles, 'NumColumns', 3, 'Box', 'off', 'FontSize', 8);
    Lgd.Layout.Tile = 'south';
end

function Fig = PlotAlphaWindow(AllRuns, Methods)
    Colours = FigureColours();
    Alphas = uniquetol(AllRuns.Alpha, 1e-9);
    Fig = NewPaperFigure(17.8, 10.5);
    Layout = tiledlayout(Fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for Tolerance = [1e-3, 1e-5]
        for Theta = [-1, +1]
            Ax = nexttile(Layout);
            hold(Ax, 'on');
            for iMethod = 1:height(Methods)
                Rows = SweepRows(AllRuns, Theta, Methods.SweepName(iMethod), Tolerance);
                for iAlpha = 1:numel(Alphas)
                    Run = RunAtAlpha(Rows, Alphas(iAlpha));
                    if isempty(Run)
                        plot(Ax, iAlpha, iMethod, '.', 'Color', Colours.NotRunDot, 'MarkerSize', 6);
                    elseif Run.Failed
                        plot(Ax, iAlpha, iMethod, 'x', 'Color', Colours.Ink, 'MarkerSize', 7, 'LineWidth', 1.3);
                    elseif Run.MaxPenetration_h <= Tolerance
                        plot(Ax, iAlpha, iMethod, 'o', 'MarkerFaceColor', Colours.Ink, ...
                             'MarkerEdgeColor', Colours.Ink, 'MarkerSize', 6);
                    else
                        plot(Ax, iAlpha, iMethod, 'o', 'MarkerFaceColor', [1, 1, 1], ...
                             'MarkerEdgeColor', Colours.Ink, 'MarkerSize', 6, 'LineWidth', 1);
                    end
                end
            end
            set(Ax, 'YDir', 'reverse', 'XLim', [0.4, numel(Alphas) + 0.6], 'YLim', [0.4, height(Methods) + 0.6], ...
                'XTick', 1:numel(Alphas), 'XTickLabel', AlphaTickLabels(Alphas), ...
                'YTick', 1:height(Methods), 'YTickLabel', cellstr(Methods.DisplayName), 'TickLength', [0, 0]);
            StyleAxes(Ax, Colours);
            if Theta == +1
                Ax.YTickLabel = {};
            end
            if Tolerance < 1e-3
                xlabel(Ax, 'penalty parameter $\alpha$', 'Color', Colours.Ink);
            end
            title(Ax, "$\theta=" + sprintf('%+d', Theta) + "$, tolerance $\varepsilon=10^{" + string(round(log10(Tolerance))) + "}$", ...
                  'FontWeight', 'normal', 'FontSize', 8.5, 'Color', Colours.Ink);
            Ax.TitleHorizontalAlignment = 'left';
            if Tolerance == 1e-3 && Theta == -1
                LegendAxes = Ax;
            end
        end
    end
    LegendHandles = [ ...
        plot(LegendAxes, NaN, NaN, 'o', 'MarkerFaceColor', Colours.Ink, 'MarkerEdgeColor', Colours.Ink, ...
             'MarkerSize', 6, 'DisplayName', 'max penetration $\leq \varepsilon h$'), ...
        plot(LegendAxes, NaN, NaN, 'o', 'MarkerFaceColor', [1, 1, 1], 'MarkerEdgeColor', Colours.Ink, ...
             'MarkerSize', 6, 'DisplayName', 'converged, penetration $> \varepsilon h$'), ...
        plot(LegendAxes, NaN, NaN, 'x', 'Color', Colours.Ink, 'MarkerSize', 7, 'LineWidth', 1.3, ...
             'DisplayName', 'not converged'), ...
        plot(LegendAxes, NaN, NaN, '.', 'Color', Colours.NotRunDot, 'MarkerSize', 6, 'DisplayName', 'not run')];
    Lgd = legend(LegendAxes, LegendHandles, 'NumColumns', 4, 'Box', 'off', 'FontSize', 8);
    Lgd.Layout.Tile = 'south';
end

function Fig = PlotSweepHeatmap(AllRuns, Methods, Quantity)
    % method x alpha grid; cell colour = value class, cell text = value
    Colours = FigureColours();
    switch Quantity
        case "penetration"
            Column       = 'MaxPenetration_h';
            Classes      = PenetrationClasses();
            FormatValue  = @FormatPenetration;
            LegendTitle  = 'maximum penetration $/\,h$';
            CellFontSize = 6.5;
        case "iterations"
            Column       = 'NewtonIterations';
            Classes      = IterationClasses();
            FormatValue  = @(Value) sprintf('%d', round(Value));
            LegendTitle  = 'total Newton iterations over 40 load steps';
            CellFontSize = 7.5;
    end
    Curves = CurveList(Methods);
    Alphas = uniquetol(AllRuns.Alpha, 1e-9);
    Fig = NewPaperFigure(17.8, 13.5);
    Layout = tiledlayout(Fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for Theta = [-1, +1]
        Ax = nexttile(Layout);
        hold(Ax, 'on');
        for iRow = 1:numel(Curves)
            Rows = SweepRows(AllRuns, Theta, Curves(iRow).SweepName, Curves(iRow).PenetrationTol);
            for iAlpha = 1:numel(Alphas)
                Run = RunAtAlpha(Rows, Alphas(iAlpha));
                if isempty(Run)
                    continue
                end
                if Run.Failed
                    Fill = Colours.NotConvergedFill;
                    CellText = 'n.c.';
                    TextColour = Colours.InkSecondary;
                else
                    Value = Run.(Column);
                    Fill = ClassFill(Value, Classes);
                    CellText = FormatValue(Value);
                    TextColour = TextColourOn(Fill, Colours.Ink);
                end
                rectangle(Ax, 'Position', [iAlpha - 0.47, iRow - 0.46, 0.94, 0.92], 'FaceColor', Fill, 'EdgeColor', 'none');
                text(Ax, iAlpha, iRow, CellText, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', CellFontSize, 'Color', TextColour);
            end
        end
        if Theta == +1
            text(Ax, 2.5, (numel(Curves) + 1) / 2, 'not run ($\alpha<1$)', 'HorizontalAlignment', 'center', ...
                 'Color', Colours.InkSecondary, 'FontSize', 8);
        end
        set(Ax, 'YDir', 'reverse', 'XLim', [0.5, numel(Alphas) + 0.5], 'YLim', [0.5, numel(Curves) + 0.5], ...
            'XTick', 1:numel(Alphas), 'XTickLabel', AlphaTickLabels(Alphas), ...
            'YTick', 1:numel(Curves), 'YTickLabel', cellstr([Curves.DisplayName]), 'TickLength', [0, 0]);
        StyleAxes(Ax, Colours);
        title(Ax, PanelTitle(Theta), 'FontWeight', 'normal', 'FontSize', 9, 'Color', Colours.Ink);
        Ax.TitleHorizontalAlignment = 'left';
        if Theta == -1
            LegendAxes = Ax;
        else
            xlabel(Ax, 'penalty parameter $\alpha$', 'Color', Colours.Ink);
        end
    end
    LegendHandles = gobjects(0);
    for iClass = 1:numel(Classes)
        LegendHandles(end + 1) = patch(LegendAxes, NaN, NaN, Classes(iClass).Fill, 'EdgeColor', 'none', ...
                                       'DisplayName', Classes(iClass).Label); %#ok<AGROW>
    end
    LegendHandles(end + 1) = patch(LegendAxes, NaN, NaN, Colours.NotConvergedFill, 'EdgeColor', 'none', ...
                                   'DisplayName', 'not converged (n.c.)');
    Lgd = legend(LegendAxes, LegendHandles, 'NumColumns', numel(LegendHandles), 'Box', 'off', 'FontSize', 8);
    Lgd.Title.String = LegendTitle;
    Lgd.Title.FontWeight = 'normal';
    Lgd.Layout.Tile = 'south';
end


%% ======================================================================
%  LaTeX tables
%  ======================================================================
function WriteSweepLatexTables(AllRuns, Methods, OutputFile)
    % needs \usepackage{booktabs} and \usepackage[table]{xcolor}; cell colours can be switched off with
    % \renewcommand{\SweepCell}[3]{#3} after \input of this file
    Curves = CurveList(Methods);
    Lines = strings(0, 1);
    Lines(end + 1, 1) = "% generated by SweepFigures.m from Sweep/data/*.mat";
    Lines(end + 1, 1) = "% needs \usepackage{booktabs} and \usepackage[table]{xcolor}";
    for Classes = {PenetrationClasses(), IterationClasses()}
        for Class = Classes{1}
            Lines(end + 1, 1) = "\definecolor{" + Class.ColourName + "}{HTML}{" + extractAfter(upper(Class.Hex), 1) + "}"; %#ok<AGROW>
        end
    end
    Lines(end + 1, 1) = "\definecolor{SweepNotConverged}{HTML}{E4E3DF}";
    Lines(end + 1, 1) = "\providecommand{\SweepCell}[3]{\cellcolor{#1}\textcolor{#2}{#3}}";
    Lines(end + 1, 1) = "";
    Lines = [Lines; SummaryTableLines(AllRuns, Methods)];
    Lines(end + 1, 1) = "";
    Lines = [Lines; DataTableLines(AllRuns, Curves, "penetration")];
    Lines(end + 1, 1) = "";
    Lines = [Lines; DataTableLines(AllRuns, Curves, "iterations")];
    FileId = fopen(OutputFile, 'w');
    fprintf(FileId, '%s\n', Lines);
    fclose(FileId);
end

function Lines = SummaryTableLines(AllRuns, Methods)
    % widest contiguous alpha block with convergence and max penetration <= tolerance, and iteration ranges
    Lines = strings(0, 1);
    Lines(end + 1, 1) = "\begin{table}[htbp]";
    Lines(end + 1, 1) = "\centering";
    Lines(end + 1, 1) = "\caption{Widest contiguous range of $\alpha$ in which the run converges with maximum penetration $\le \varepsilon h$, and total Newton iterations over all converged runs. Tested ranges: $\alpha\in[10^{-2},10^{3}]$ for $\theta=-1$ and $\alpha\in[1,10^{3}]$ for $\theta=+1$; Penalty and augmented Lagrange do not depend on $\theta$.}";
    Lines(end + 1, 1) = "\label{tab: sweep summary}";
    Lines(end + 1, 1) = "\small";
    Lines(end + 1, 1) = "\begin{tabular}{l ccc ccc}";
    Lines(end + 1, 1) = "\toprule";
    Lines(end + 1, 1) = " & \multicolumn{3}{c}{$\theta=-1$} & \multicolumn{3}{c}{$\theta=+1$} \\";
    Lines(end + 1, 1) = "\cmidrule(lr){2-4}\cmidrule(lr){5-7}";
    Lines(end + 1, 1) = "Formulation & $\varepsilon=10^{-3}$ & $\varepsilon=10^{-5}$ & Iterations & $\varepsilon=10^{-3}$ & $\varepsilon=10^{-5}$ & Iterations \\";
    Lines(end + 1, 1) = "\midrule";
    HasFurtherBlocks = false;
    for iMethod = 1:height(Methods)
        Row = Methods.DisplayName(iMethod);
        for Theta = [-1, +1]
            Iterations = [];
            for Tolerance = [1e-3, 1e-5]
                Rows = SweepRows(AllRuns, Theta, Methods.SweepName(iMethod), Tolerance);
                IsCorrect = ~Rows.Failed & Rows.MaxPenetration_h <= Tolerance;
                [Text, IsWidestOnly] = WidestAlphaBlock(Rows.Alpha, IsCorrect);
                HasFurtherBlocks = HasFurtherBlocks || ~IsWidestOnly;
                Row = Row + " & " + Text;
                Iterations = [Iterations; Rows.NewtonIterations(~Rows.Failed)]; %#ok<AGROW>
            end
            Row = Row + " & " + sprintf('%d--%d', min(Iterations), max(Iterations));
        end
        Lines(end + 1, 1) = Row + " \\"; %#ok<AGROW>
    end
    Lines(end + 1, 1) = "\bottomrule";
    Lines(end + 1, 1) = "\end{tabular}";
    if HasFurtherBlocks
        Lines(end + 1, 1) = "";
        Lines(end + 1, 1) = "\smallskip";
        Lines(end + 1, 1) = "\footnotesize $^{\dagger}$ the tolerance is also met at isolated values of $\alpha$ outside this range.";
    end
    Lines(end + 1, 1) = "\end{table}";
end

function [Text, IsWidestOnly] = WidestAlphaBlock(Alpha, IsOk)
    % same block search as AlphaRanges in SweepAlpha_call.m
    Change     = diff([0; IsOk(:); 0]);
    BlockStart = find(Change == 1);
    BlockEnd   = find(Change == -1) - 1;
    IsWidestOnly = numel(BlockStart) <= 1;
    if isempty(BlockStart)
        Text = "--";
        return
    end
    [~, Widest] = max(log10(Alpha(BlockEnd) ./ Alpha(BlockStart)));
    Text = "$[" + AlphaPower(Alpha(BlockStart(Widest))) + "," + AlphaPower(Alpha(BlockEnd(Widest))) + "]";
    if ~IsWidestOnly
        Text = Text + "^{\dagger}";
    end
    Text = Text + "$";
end

function Lines = DataTableLines(AllRuns, Curves, Quantity)
    % rows = alpha (grouped by theta), columns = formulations; cells coloured like the heatmaps
    switch Quantity
        case "penetration"
            Column  = 'MaxPenetration_h';
            Classes = PenetrationClasses();
            Caption = "Maximum penetration $/\,h$ at the end of loading for every run of the sweep. Cell colours: four classes bounded by $10^{-8}$ and the tolerances $\varepsilon=10^{-5}$, $10^{-3}$; n.c.\ = not converged.";
            Label   = "tab: sweep penetration";
        case "iterations"
            Column  = 'NewtonIterations';
            Classes = IterationClasses();
            Caption = "Total number of Newton iterations over the 40 load steps for every run of the sweep, including the outer-loop iterations of the augmented and penalty--Nitsche formulations; n.c.\ = not converged.";
            Label   = "tab: sweep iterations";
    end
    Alphas = uniquetol(AllRuns.Alpha, 1e-9);
    Lines = strings(0, 1);
    Lines(end + 1, 1) = "\begin{table}[htbp]";
    Lines(end + 1, 1) = "\centering";
    Lines(end + 1, 1) = "\caption{" + Caption + "}";
    Lines(end + 1, 1) = "\label{" + Label + "}";
    Lines(end + 1, 1) = "\small";
    Lines(end + 1, 1) = "\setlength{\tabcolsep}{4pt}";
    Lines(end + 1, 1) = "\begin{tabular}{c" + join(repmat("c", 1, numel(Curves)), "") + "}";
    Lines(end + 1, 1) = "\toprule";
    Lines(end + 1, 1) = "$\alpha$ & Penalty & Nitsche & Augmented & \multicolumn{2}{c}{Penalty--Nitsche} & Augmented \\";
    Lines(end + 1, 1) = "\cmidrule(lr){5-6}";
    Lines(end + 1, 1) = " & & & Nitsche & $\varepsilon=10^{-3}$ & $\varepsilon=10^{-5}$ & Lagrange \\";
    for Theta = [-1, +1]
        Lines(end + 1, 1) = "\midrule"; %#ok<AGROW>
        Lines(end + 1, 1) = "\multicolumn{" + string(numel(Curves) + 1) + "}{l}{" + extractAfter(PanelTitle(Theta), 4) + "} \\"; %#ok<AGROW>
        for iAlpha = 1:numel(Alphas)
            Row = "$" + AlphaPower(Alphas(iAlpha)) + "$";
            HasRuns = false;
            for Curve = Curves
                Run = RunAtAlpha(SweepRows(AllRuns, Theta, Curve.SweepName, Curve.PenetrationTol), Alphas(iAlpha));
                if isempty(Run)
                    Row = Row + " & ";
                    continue
                end
                HasRuns = true;
                if Run.Failed
                    Row = Row + " & \SweepCell{SweepNotConverged}{black}{n.c.}";
                else
                    Value = Run.(Column);
                    Class = Classes(find(Value <= [Classes.UpperBound], 1));
                    if Quantity == "penetration"
                        CellText = FormatPenetration(Value);
                    else
                        CellText = sprintf('%d', round(Value));
                    end
                    Row = Row + " & \SweepCell{" + Class.ColourName + "}{" + Class.TextColourName + "}{" + CellText + "}";
                end
            end
            if HasRuns
                Lines(end + 1, 1) = Row + " \\"; %#ok<AGROW>
            end
        end
    end
    Lines(end + 1, 1) = "\bottomrule";
    Lines(end + 1, 1) = "\end{tabular}";
    Lines(end + 1, 1) = "\end{table}";
end


%% ======================================================================
%  Classes, colours, formatting
%  ======================================================================
function Classes = PenetrationClasses()
    % ordinal classes of max penetration / h (blue ramp steps 250/400/550/700)
    Classes = MakeClasses([1e-8, 1e-5, 1e-3, Inf], ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"], ...
        ["$\leq 10^{-8}$ (practically exact)", "$10^{-8}$ to $10^{-5}$", "$10^{-5}$ to $10^{-3}$", "$> 10^{-3}$"], ...
        "SweepPenetration");
end

function Classes = IterationClasses()
    % ordinal classes of total Newton iterations (orange ramp built on categorical orange)
    Classes = MakeClasses([130, 160, 250, Inf], ["#f19d70", "#eb6834", "#bf4a1a", "#84300c"], ...
        ["$\leq 130$", "131 to 160", "161 to 250", "$> 250$"], "SweepIterations");
end

function Classes = MakeClasses(UpperBounds, HexColours, Labels, ColourPrefix)
    Ink = HexToRgb("#0b0b0b");
    Classes = struct('UpperBound', {}, 'Hex', {}, 'Fill', {}, 'Label', {}, 'ColourName', {}, 'TextColourName', {});
    for iClass = 1:numel(UpperBounds)
        Fill = HexToRgb(HexColours(iClass));
        if isequal(TextColourOn(Fill, Ink), Ink)
            TextColourName = "black";
        else
            TextColourName = "white";
        end
        Classes(iClass) = struct('UpperBound', UpperBounds(iClass), 'Hex', HexColours(iClass), 'Fill', Fill, ...
            'Label', Labels(iClass), 'ColourName', ColourPrefix + string(iClass), 'TextColourName', TextColourName);
    end
end

function Fill = ClassFill(Value, Classes)
    Fill = Classes(find(Value <= [Classes.UpperBound], 1)).Fill;
end

function Colours = FigureColours()
    Colours.Ink              = HexToRgb("#0b0b0b");
    Colours.InkSecondary     = HexToRgb("#52514e");
    Colours.Grid             = HexToRgb("#e4e3df");
    Colours.NotRunFill       = HexToRgb("#f1f0ec");
    Colours.NotRunDot        = HexToRgb("#c9c8c2");
    Colours.NotConvergedFill = HexToRgb("#e4e3df");
end

function Rgb = HexToRgb(Hex)
    Hex = char(Hex);
    Rgb = reshape(sscanf(Hex(2:end), '%2x'), 1, 3) / 255;
end

function TextColour = TextColourOn(Fill, Ink)
    % ink or white, whichever has the higher WCAG contrast on the fill
    ContrastInk   = (RelativeLuminance(Fill) + 0.05) / (RelativeLuminance(Ink) + 0.05);
    ContrastWhite = 1.05 / (RelativeLuminance(Fill) + 0.05);
    if ContrastInk >= ContrastWhite
        TextColour = Ink;
    else
        TextColour = [1, 1, 1];
    end
end

function Luminance = RelativeLuminance(Rgb)
    Linear = Rgb / 12.92;
    IsHigh = Rgb > 0.03928;
    Linear(IsHigh) = ((Rgb(IsHigh) + 0.055) / 1.055) .^ 2.4;
    Luminance = 0.2126 * Linear(1) + 0.7152 * Linear(2) + 0.0722 * Linear(3);
end

function Label = FormatPenetration(Value)
    % 3.5e-4 -> '$3.5\cdot10^{-4}$'; one more digit only if rounding would print a class boundary
    % for a value lying on its other side (e.g. 1.03e-5 must not print as 1.0e-5)
    if Value == 0
        Label = '$0$';
        return
    end
    Boundaries = [1e-8, 1e-5, 1e-3];
    Exponent = floor(log10(Value));
    Digits = 1;
    Rounded = round(Value / 10^Exponent, Digits) * 10^Exponent;
    if any(abs(Rounded - Boundaries) <= 1e-12 * Boundaries & Value ~= Boundaries)
        Digits = 2;
    end
    Mantissa = round(Value / 10^Exponent, Digits);
    if Mantissa >= 10
        Mantissa = Mantissa / 10;
        Exponent = Exponent + 1;
    end
    Label = sprintf('$%.*f\\cdot10^{%d}$', Digits, Mantissa, Exponent);
end

function Text = AlphaPower(Alpha)
    % 0.0316 -> '10^{-1.5}'
    Text = "10^{" + string(round(2 * log10(Alpha)) / 2) + "}";
end

function Labels = AlphaTickLabels(Alphas)
    Powers = round(2 * log10(Alphas(:))) / 2;
    Labels = cellstr("$10^{" + string(Powers) + "}$");
end

function Title = PanelTitle(Theta)
    if Theta == -1
        Title = "(a) $\theta=-1$ (skew-symmetric)";
    else
        Title = "(b) $\theta=+1$ (symmetric)";
    end
end

function Shift = FailedMarkerShift(SweepName)
    % Penalty and Augmented Lagrange both fail at alpha = 1e3; move their "not converged" markers apart
    Shift = 1;
    if SweepName == "Penalty"
        Shift = 0.84;
    elseif SweepName == "Augumented Lagrange"
        Shift = 1.19;
    end
end

function Fig = NewPaperFigure(WidthCm, HeightCm)
    Fig = figure('Units', 'centimeters', 'Position', [2, 2, WidthCm, HeightCm], 'Color', 'w', ...
        'DefaultTextInterpreter', 'latex', 'DefaultAxesTickLabelInterpreter', 'latex', ...
        'DefaultLegendInterpreter', 'latex', 'DefaultAxesFontSize', 9, 'DefaultTextFontSize', 9);
end

function StyleAxes(Ax, Colours)
    set(Ax, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.6, 'FontSize', 9, ...
        'XColor', Colours.InkSecondary, 'YColor', Colours.InkSecondary, ...
        'GridColor', Colours.Grid, 'GridAlpha', 1, 'Layer', 'bottom');
end