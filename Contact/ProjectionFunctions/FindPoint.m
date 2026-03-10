function Result = FindPoint(Body,Point)

    DofsAtNode = Body.DofsAtNode;
    SurfaceNodes = Body.contact.nodalid;
    nloc = Body.nloc;
     
    % current position of the "possible contact" nodes of the Body
    SurfacePoints_X = Body.q(xlocChosen(DofsAtNode,SurfaceNodes,1)) + ... % coords on X axis;
                      Body.u(xlocChosen(DofsAtNode,SurfaceNodes,1));
    
    SurfacePoints_Y = Body.q(xlocChosen(DofsAtNode,SurfaceNodes,2)) + ... % coords on Y axis
                      Body.u(xlocChosen(DofsAtNode,SurfaceNodes,2));
    
    SurfacePoints = [SurfacePoints_X SurfacePoints_Y]; % nodes of the contact surfaces of the contact body 
    
    distances = vecnorm(SurfacePoints - Point, 2, 2); % distances between all target contact node and the point
        
    [~, sortedIndices] = sort(distances); % sorting and choosing two closest
    a = SurfaceNodes(sortedIndices(1)); 
    b = SurfaceNodes(sortedIndices(2)); 

    position_a = SurfacePoints(sortedIndices(1),:); 
    position_b = SurfacePoints(sortedIndices(2),:); 

    [xy, distance, t_a]  = distance2curve([position_a; position_b], Point, 'linear'); 
    tol = 1e-6; % Tolerance for error margin 
   
    info = []; % array, where we will store info of the contact point projection 

    if ~( (t_a < tol || abs(t_a - 1) < tol) && distance > tol )  % sanity check that the point isn't outside  (t_a ~= 0, 1)
                                                                 % and having distance > 0 at the same time                
        % To what element these nodes belong
        % this is a linear element, therefore two nodes are uniquely belong to one element only 
        ElemenNumber = find(any(nloc == a, 2) & any(nloc == b, 2)); 

        if ~isempty(ElemenNumber) % sanity check (probably need to be removed)
            [X,U] = GetCoorDisp(ElemenNumber,nloc,Body.P0,Body.u); % position of nodes of the element                 
            central = Nm_2412(0,0)*(X+U); % Position of the central point of the chosen element     
            % Finding external normal to the element (a central point helps identify the outward direction)
            Normal = Normal3points(central,position_a',position_b'); % outwards normal of element (algorithm doesn't depend on the order a and b)            
            info = [info; xy Normal ElemenNumber];
        end             
        
    end 

    Result.Gap = 0; % we always have something to work with;
    
    if isempty(info) == false % we actually have contact
        
        % info extraction
        Position = info(:,1:2);
        Normal = info(:,3:4)';

        % Saving
        Result.Index = info(:,5); % element on which point projected       
        Result.Gap = (Point - Position) * Normal;  % gap calculation
        Result.Normal = Normal; % normal on the surf.
        Result.Position = Position'; % point projection

   end

    