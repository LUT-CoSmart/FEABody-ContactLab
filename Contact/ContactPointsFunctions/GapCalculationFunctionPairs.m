function Gap = GapCalculationFunctionPairs(ContactBody,TargetBody,ContactPointsName,n)
    
    ContactPoints = ContactPointsFunction(ContactBody, ContactPointsName, n);
    m = size(ContactPoints,1);
    Gap = zeros(m,1); 
    %% TODO: make it over all points simultaneously, working with array or in parallel  
    for ii = 1:m % loop over all contact points
      
        ContactPoint = ContactPoints(ii,:);
     
        % Search the attributes of the corresponding point on the target surface 
        Outcome = FindPoint(TargetBody,ContactPoint);     
        Gap(ii) = Outcome.Gap;
    end