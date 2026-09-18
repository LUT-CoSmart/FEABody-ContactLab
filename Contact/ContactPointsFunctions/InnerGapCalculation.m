function Gap = InnerGapCalculation(ContactBody,TargetBody,ContactPointsName,n)
    
    Gap = 0;    
    ContactPoints = ContactPointsFunction(ContactBody, ContactPointsName, n);
    
    %% TODO: make it over all points simultaneously, working with array or in parallel  
    %% possible if yo vectorize "FindPoint" 
    for ii = 1:size(ContactPoints,1) % loop over all contact points
      
        ContactPoint = ContactPoints(ii,:);
     
        % Search the attributes of the corresponding point on the target surface 
        Outcome = FindPoint(TargetBody,ContactPoint);
        if Outcome.Gap < 0 
            Gap = Gap + abs(Outcome.Gap);
        end
    end