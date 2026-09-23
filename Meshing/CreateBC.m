
function Body = CreateBC(Body,type)

    DofsAtNode = Body.DofsAtNode;
    
     % Set vector of linear constraints
    if ~isfield(Body,'bc') || isempty(Body.bc)
        bc = true(1,Body.nx); % Initialize BC only on the first call      
    else
        bc = Body.bc;   % keep prescribed BCs if any 
    end

    if ~strcmp(type, 'none') 
        NodalIDAt = FindGlobNodalID(Body.P0,Body.loc,Body.shift); % let's find global nodal ID's x=0
    end
    
    if strcmp(type, 'all')
        bcInd=xlocChosen(DofsAtNode,NodalIDAt,1:2);  % all fixed at clambed end
    elseif strcmp(type, 'ux')
        bcInd=xlocChosen(DofsAtNode,NodalIDAt,1);
    elseif strcmp(type, 'uy')
        bcInd=xlocChosen(DofsAtNode,NodalIDAt,2);
    elseif strcmp(type, 'none')
        bcInd=[];
    else
        error('Chose the right type of BC!')
    end    

    bcInd = bcInd(bcInd ~= 0); % Add new constraints
    
    if bcInd~=0
        bc(bcInd)=0;                       % number of degrees of freedom of system after linear constraints  
    end
    
    Body.bc = bc;
    Body.ndof = sum(bc);     % Number of unconstrained DOFs  