function [ContactPointfunc, Gapfunc, GapfuncPairs] = ContactPointSetting(PointsofInterest,approachBasis)
        
    ContactPointsName = PointsofInterest.Name;
    n = PointsofInterest.n;
         
    if ContactPointsName == "Gauss"
        [z,w] = gauleg2(-1,1,n);
    elseif ContactPointsName == "LinSpace"
        n = max(n,2); % sanity check goes here
        z = geospace(-1,1,n)';
        w = diff(z);
        z = z(2:end);
    else     
        if approachBasis == "Lagrange"
           z = 1;
           w = 2;
        else
           z = [-1; 1];
           w = [ 1; 1];
        end   
    end

    ContactPointfunc =  @(ContactBody) ContactPointsFunction(ContactBody,z,w);
    Gapfunc = @(ContactBody,TargetBody) InnerGapCalculation(ContactBody,TargetBody,z,w); 
    GapfuncPairs = @(ContactBody,TargetBody) GapCalculationFunctionPairs(ContactBody,TargetBody,z,w);