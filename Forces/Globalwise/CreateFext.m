function Body = CreateFext(currentStep,Nsteps,Body,type)

    DofsAtNode = Body.DofsAtNode;
    Fext = zeros(Body.nx,1); % Initialize vector of ext forces
    
    % Force check
    if ~isfield(Body, 'Fext')
        Body.Fext.x = 0;
        Body.Fext.y = 0;
        Body.Fext.loc.x = 'all';        
        Body.Fext.loc.y = 'all';
    else
        if ~isfield(Body.Fext.loc, 'x')
           Body.Fext.loc.x = 'all';
           Body.Fext.x = 0;
        end
        if ~isfield(Body.Fext.loc, 'y')
           Body.Fext.loc.y = 'all';
           Body.Fext.y = 0;
        end
    end

    NodalIDAt = FindGlobNodalID(Body.P0,Body.Fext.loc,Body.shift);
    NumberOfIDs = length(NodalIDAt);

    switch type
           case "linear"
                Loadstep =  currentStep/Nsteps;
           case "exponential"     
                coef = 5;
                Loadstep = (exp((currentStep/Nsteps).^coef) - 1) / (exp(1) - 1);
           case "quadratic"
                Loadstep = (currentStep/Nsteps)^2;
           case "cubic"
                Loadstep = (currentStep/Nsteps)^3;    
           case "quartic"  
                Loadstep = (currentStep/Nsteps)^4;    
           otherwise
                error('Unknown loading type')                 
     end  

    ByNodes = max(NumberOfIDs - 1,1);
    
    Load_x = Body.Fext.x*Loadstep/ByNodes;
    Load_y = Body.Fext.y*Loadstep/ByNodes;
    
    Fext(xlocChosen(DofsAtNode,NodalIDAt,1)) = Load_x;
    Fext(xlocChosen(DofsAtNode,NodalIDAt,2)) = Load_y;
    
    if NumberOfIDs > 1
    
        % Endpoints of the ordered loaded edge.
        EdgeNodes = NodalIDAt([1,end]);
    
        Fext(xlocChosen(DofsAtNode,EdgeNodes,1)) = Load_x/2;
        Fext(xlocChosen(DofsAtNode,EdgeNodes,2)) = Load_y/2;
    
    end

    Body.Fext.vec = Fext;
    Body.Fext.nodalid = NodalIDAt;
