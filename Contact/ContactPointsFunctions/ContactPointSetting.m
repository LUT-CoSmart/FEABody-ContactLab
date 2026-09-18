function [ContactPointfunc, Gapfunc, GapfuncPairs] = ContactPointSetting(PointsofInterest)
    
    addpath("Contact\ProjectionFunctions")   
    ContactPoints = PointsofInterest.Name;
    n = PointsofInterest.n;
  
    % Sanity check
    if  ContactPoints == "LinSpace"
            if n < 2 
               warning("Number of points is not enough, it is set to 2 ")
               n = 2; 
            end     
    end
    
    ContactPointfunc =  @(ContactBody) ContactPointsFunction(ContactBody,ContactPoints,n);
    Gapfunc = @(ContactBody,TargetBody) InnerGapCalculation(ContactBody,TargetBody,ContactPoints,n); 
    GapfuncPairs = @(ContactBody,TargetBody) GapCalculationFunctionPairs(ContactBody,TargetBody,ContactPoints,n);