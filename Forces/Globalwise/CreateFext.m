
function Body = CreateFext(currentStep,Nsteps,Body,type)
    
    DofsAtNode = Body.DofsAtNode;
    Fext = zeros(Body.nx,1); % Initialize vector of ext forces
    
    NodalIDAt = FindGlobNodalID(Body.P0,Body.Fext.loc,Body.shift);
    NumberOfIDs = length(NodalIDAt);
    ByNodes =  NumberOfIDs-1;

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
     
    
    Load_x = Body.Fext.x*Loadstep/ByNodes;
    Load_y = Body.Fext.y*Loadstep/ByNodes;
    
    
    Fext(xlocChosen(DofsAtNode,NodalIDAt,1))= Load_x;       
    Fext(xlocChosen(DofsAtNode,NodalIDAt,2))= Load_y;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NodalIDAtEdge1 = FindGlobNodalID(Body.P0,Body.edge1.loc,Body.shift);
    NodalIDAtEdge2 = FindGlobNodalID(Body.P0,Body.edge2.loc,Body.shift);

    nodeEdge1 = intersect(NodalIDAt, NodalIDAtEdge1);
    nodeEdge2 = intersect(NodalIDAt, NodalIDAtEdge2);

    Fext(xlocChosen(DofsAtNode,nodeEdge1,1))= Load_x/2;
    Fext(xlocChosen(DofsAtNode,nodeEdge1,2))= Load_y/2;

    Fext(xlocChosen(DofsAtNode,nodeEdge2,1))= Load_x/2;
    Fext(xlocChosen(DofsAtNode,nodeEdge2,2))= Load_y/2;

    
    Body.Fext.vec = Fext;
    Body.Fext.nodalid = NodalIDAt;