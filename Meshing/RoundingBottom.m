function Body = RoundingBottom(Body, Radius, ArcCenterX)

    % Lifts the nodes so that the bottom edge follows a circular arc; x is not changing
    % ArcCenterX (optional): x of the lowest point of the arc; default is the middle of the body
    % Idea is to consider sym of the Hertz prblem, to save computational time
    if nargin < 3
        ArcCenterX = Body.shift.x + Body.Lx/2;
    end
    LargestDistance = max(abs(Body.P0(:,1) - ArcCenterX));
    if Radius < LargestDistance
        warning('Radius must be at least the largest horizontal distance from the arc centre, set to it');
        Radius = LargestDistance;
    end

    BottomX = ArcCenterX;
    BottomY = Body.shift.y;

    x = Body.P0(:,1);
    y = Body.P0(:,2);
    ArcRise = Radius - sqrt(Radius^2 - (x - BottomX).^2);
    ArcRise = ArcRise - min(ArcRise);      
    
    Body.P0(:,2) = y + ArcRise .* (1 - (y - BottomY)/Body.Ly);
    Body.q = reshape(Body.P0.', [], 1);