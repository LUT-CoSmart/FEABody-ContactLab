function Result = FindPoint(Body,Point)
    DofsAtNode = Body.DofsAtNode;
    SurfaceNodes = Body.contact.nodalid;
    nloc = Body.nloc;
    % current position of the surface nodes
    idX = xlocChosen(DofsAtNode,SurfaceNodes,1);
    idY = xlocChosen(DofsAtNode,SurfaceNodes,2);
    SurfacePoints_X = Body.q(idX) + Body.u(idX); % coords on X axis;
    SurfacePoints_Y = Body.q(idY) + Body.u(idY); % coords on Y axis
    SurfacePoints = [SurfacePoints_X SurfacePoints_Y]; % nodes of the contact surfaces of the contact body
    distances = vecnorm(SurfacePoints - Point, 2, 2); % distances between all target contact node and the point
    [~, sortedIndices] = sort(distances); % sorting and choosing two closest
    a = SurfaceNodes(sortedIndices(1));
    b = SurfaceNodes(sortedIndices(2));
    position_a = SurfacePoints(sortedIndices(1),:);
    position_b = SurfacePoints(sortedIndices(2),:);
    tol = 1e-6; % Tolerance for error margin
    [xy, ~, ~, t_raw] = ProjectToSegment(position_a, position_b, Point);

    % projection is outside [a,b] -> try the segment on the other side of a
    if t_raw < -tol || t_raw > 1 + tol
        ia = sortedIndices(1);
        ib = ia - sign(sortedIndices(2) - ia);
        if ib >= 1 && ib <= numel(SurfaceNodes)
            b = SurfaceNodes(ib);
            position_b = SurfacePoints(ib,:);
            [xy, ~, ~, t_raw] = ProjectToSegment(position_a, position_b, Point);
        end
    end

    info = []; 
    if t_raw >= -tol && t_raw <= 1 + tol   % node vertices (t = 0 or 1) are now accepted
        ElemenNumber = find(any(nloc == a, 2) & any(nloc == b, 2));
        if ~isempty(ElemenNumber) % sanity check (probably need to be removed)
            [X,U] = GetCoorDisp(ElemenNumber,nloc,Body.P0,Body.u); % position of nodes of the element
            central = Nm_2412(0,0)*(X+U); % Position of the central point of the chosen element
            Normal = Normal3points(central,position_a',position_b'); % outwards normal of element
            info = [info; xy Normal ElemenNumber];
        end
    end
    Result.Gap = 0; % we always have something to work with;
    if isempty(info) == false % we actually have contact
        Position = info(:,1:2);
        Normal = info(:,3:4)';
        Result.Index = info(:,5); % element on which point projected
        Result.Gap = (Point - Position) * Normal;  % gap calculation
        Result.Normal = Normal; % normal on the surf.
        Result.Position = Position'; % point projection
    end