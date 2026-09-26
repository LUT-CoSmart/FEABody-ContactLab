function [count,ContactGeometry, TargetGeometry, Gap, Normal] = Projection(ContactPoints,ContactPointsElements,ContactBody,TargetBody)

    nloc_cont = ContactBody.nloc;
    xloc_cont = ContactBody.xloc;

    count = 0;

    Gap = [];
    Normal = [];

    ContactDofs = [];
    ContactCoordsXi = [];
    ContactCoords = [];
    ContactDisp = [];
    ContactIndex = []; 
    
    TargetDofs = [];
    TargetCoordsXi = [];
    TargetCoords = [];
    TargetDisp = [];
       
    %% TODO: make it over all points simultaneously, working with array or in parallel 
    for ii = 1:size(ContactPoints,1) % loop over all contact points
                
        ContactPoint = ContactPoints(ii,:);  
        Outcome = FindPoint(TargetBody,ContactPoint); % forward projection
        nloc_targ = TargetBody.nloc;
        
        % Checking the condition of the penalty approach
        if Outcome.Gap < 0 % we have meaningful results after the search         
            % negative because after penetration direction is from contact point
            % to its projection, which is an opposite to the outwards normal
            
            count = count + 1; % number of contact points
            gap = abs(Outcome.Gap);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ElemenNumber = ContactPointsElements(ii);                
            DOFpositions_cont = xloc_cont(ElemenNumber,:); % DOFs of the contact element under consideration
            [X_cont,U_cont] = GetCoorDisp(ElemenNumber,nloc_cont,ContactBody.P0,ContactBody.u);
            [xi_cont,eta_cont] = FindIsoCoord(X_cont,U_cont,ContactPoint');  % finding the isocoordinates of the contact point
            vec_cont = [xi_cont; eta_cont];
            
            % index - element under consideration
            Index = Outcome.Index;
            DOFpositions_targ = TargetBody.xloc(Index,:);
            [X_targ, U_targ] = GetCoorDisp(Index,nloc_targ,TargetBody.P0,TargetBody.u);
            [xi_targ,eta_targ] = FindIsoCoord(X_targ,U_targ,Outcome.Position); % finding isoparametric coodinates of the target point    
            vec_targ = [xi_targ; eta_targ];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            Gap = [Gap gap];
            Normal = [Normal Outcome.Normal]; % outwards normal (from the target body to the contact one)

            ContactDofs = [ContactDofs DOFpositions_cont(:)];
            ContactCoordsXi = [ContactCoordsXi vec_cont];
            ContactCoords = [ContactCoords X_cont];
            ContactDisp = [ContactDisp U_cont];
            ContactIndex = [ContactIndex ii]; 

            TargetDofs = [TargetDofs DOFpositions_targ(:)];
            TargetCoordsXi = [TargetCoordsXi vec_targ];
            TargetCoords = [TargetCoords X_targ];
            TargetDisp = [TargetDisp U_targ];

        end
    end
    
    ContactGeometry.Dofs = ContactDofs;
    ContactGeometry.CoordsXi = ContactCoordsXi;
    ContactGeometry.Coords = ContactCoords;    
    ContactGeometry.Disp= ContactDisp;    
    ContactGeometry.Index = ContactIndex; 

    TargetGeometry.Dofs = TargetDofs;
    TargetGeometry.CoordsXi = TargetCoordsXi; 
    TargetGeometry.Coords = TargetCoords;    
    TargetGeometry.Disp= TargetDisp;  