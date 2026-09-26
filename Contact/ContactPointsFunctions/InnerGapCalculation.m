function [Gap,GapMax] = InnerGapCalculation(ContactBody,TargetBody,ContactPointsName,n) 
    

    Gap = 0;    
    [ContactPoints,~,ContactAreas] = ContactPointsFunction(ContactBody, ContactPointsName, n);  
    
    AreaTot = 0;                                                                          
    GapMax = 0;                                                                           
    
    %% TODO: make it over all points simultaneously, working with array or in parallel  
    %% possible if yo vectorize "FindPoint" 
    for ii = 1:size(ContactPoints,1) % loop over all contact points
      
        ContactPoint = ContactPoints(ii,:);
     
        % Search the attributes of the corresponding point on the target surface 
        Outcome = FindPoint(TargetBody,ContactPoint);
        if Outcome.Gap < 0 
            Gap = Gap + abs(Outcome.Gap)*ContactAreas(ii);                               
            AreaTot = AreaTot + ContactAreas(ii);                                         
            GapMax = max(GapMax, abs(Outcome.Gap));                                       
        end
    end

    if AreaTot > 0                                                                      
        Gap = Gap/AreaTot;   % area-weighted mean penetration                    
    end                                                                                  