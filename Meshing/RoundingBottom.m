function Body = RoundingBottom(Body, Radius)

    % Lifts the nodes so that the bottom edge follows a circular arc; x is not changing
    if Radius < Body.Lx/2
        warning('Radius must be at least half of the width, set to the half');
        Radius =  Body.Lx/2;
    end    

    BottomX = Body.shift.x + Body.Lx/2;
    BottomY = Body.shift.y;

    x = Body.P0(:,1);
    y = Body.P0(:,2);
    ArcRise = Radius - sqrt(Radius^2 - (x - BottomX).^2);
    ArcRise = ArcRise - min(ArcRise);      
    
    Body.P0(:,2) = y + ArcRise .* (1 - (y - BottomY)/Body.Ly);
    Body.q = reshape(Body.P0.', [], 1);