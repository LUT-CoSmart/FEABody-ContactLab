function [xy, distance, t_a, t_raw] = ProjectToSegment(a, b, P) % SUBS FOR DISTANCE TO CURVEE - WORKS, BECAUSE WE HAVE LINEAR 
% ELEMNET HERE 
    ab = b - a;
    ap = P - a;
    L2 = sum(ab.*ab);
    if L2 == 0
        t_raw = 0;
    else
        t_raw = sum(ap.*ab) / L2;     
    end
    t_a = min(max(t_raw, 0), 1);      % same meaning as distance2curve's t_a
    xy  = a + t_a*ab;                 % closest point (1x2, like distance2curve)
    d = P - xy;
    distance = sqrt(sum(d.*d));
end