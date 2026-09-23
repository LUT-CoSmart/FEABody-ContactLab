function [ContactPoints, ContactPointsElements,ContactAreas] = ContactPointsFunction(ContactBody, z,w)
 
    % nloc_cont = ContactBody.nloc;
    DofsAtNode_cont = ContactBody.DofsAtNode;
    [ContactSegments, ContactSegmentsElements] = ContactSegmentsFunction(ContactBody);
        
    ContactPoints = []; % taking the first node position, otherwise later it will be omitted 
    ContactPointsElements = []; % taking the first node number
    ContactAreas = [];
    
    for ii = 1:size(ContactSegments,1)

        nodes = ContactSegments(ii,:);

        idX = xlocChosen(DofsAtNode_cont,nodes,1);
        idY = xlocChosen(DofsAtNode_cont,nodes,2);

        X0 = [ContactBody.q(idX), ContactBody.q(idY)];
        X = X0 + [ContactBody.u(idX), ContactBody.u(idY)];

        a = X(1,:);
        b = X(2,:);

        Points = 0.5*(1-z)*a + 0.5*(1+z)*b;

        % Reference segment length, matching your present area convention.
        L0 = norm(X0(2,:) - X0(1,:))/2;

        Areas = w*L0*ContactBody.Lz;

        ElementNumber = ContactSegmentsElements(ii);

        ContactPoints = [ContactPoints; Points];

        ContactPointsElements = [ContactPointsElements; ElementNumber*ones(numel(z),1)];

        ContactAreas = [ContactAreas; Areas];

    end

