function NodalID = FindGlobNodalID(P,loc,shift)

    tol = 2*sqrt(eps);

    % return coordinates to the square [(0.0)-(Lx,Ly)]
    x = P(:,1) - shift.x;
    y = P(:,2) - shift.y;

    % Selection along x.
    if isnumeric(loc.x) && ismember(numel(loc.x),[1,2])
        inX = x >= min(loc.x)-tol & x <= max(loc.x)+tol;
    elseif strcmp(loc.x,'all')
        inX = true(size(P,1),1);
    else
        error('loc.x must be a scalar, [a,b], or ''all''.');
    end

    % Selection along y.
    if isnumeric(loc.y) && ismember(numel(loc.y),[1,2])
        inY = y >= min(loc.y)-tol &  y <= max(loc.y)+tol;
    elseif strcmp(loc.y,'all')
        inY = true(size(P,1),1);
    else
        error('loc.y must be a scalar, [a,b], or ''all''.');
    end

    NodalID = find(inX & inY).';

end