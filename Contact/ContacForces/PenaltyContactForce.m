function Fc = PenaltyContactForce(ContactBody,TargetBody,approach)
    
    ContactPointfunc = approach.ContactPointfunc;
    penalty = approach.penalty;

    % Target Body is a body points are projected
    % Contact Body is a body points are taken for projection
    Fcont = zeros(ContactBody.nx,1);
    Ftarg = zeros(TargetBody.nx,1);

    [ContactPoints, ContactPointsElements] = ContactPointfunc(ContactBody);

    [count,ContactGeometry, TargetGeometry, Gaps, Normals] = Projection(ContactPoints,ContactPointsElements,ContactBody,TargetBody); 

    if count~=0 % we have contact
        
        for i = 1:count % loop all over contact points

            xi_cont = ContactGeometry.CoordsXi(:,i) ;    
            xi_targ = TargetGeometry.CoordsXi(:,i);
            
            Gap = Gaps(i);
            
            Normal = Normals(:,i);
            
            Normal_cont = -Normal;
            Normal_targ =  Normal; 
            
            % calculation of the forces applied to the nodes of contact elemnet 
            Fcont_loc = penalty * Gap * Normal_cont;                                                                              
            Ftarg_loc = penalty * Gap * Normal_targ; 
                       
            % Redistribution over the nodes
            DOFpositions_cont = ContactGeometry.Dofs(:,i);
            DOFpositions_targ = TargetGeometry.Dofs(:,i);

            Fcont(DOFpositions_cont) = Fcont(DOFpositions_cont) + Nm_2412(xi_cont(1),xi_cont(2))' * Fcont_loc;
            Ftarg(DOFpositions_targ) = Ftarg(DOFpositions_targ) + Nm_2412(xi_targ(1),xi_targ(2))' * Ftarg_loc;   

        end 
    end  
    
    % Assemblace
    Fc = [Fcont;Ftarg]; 

    
