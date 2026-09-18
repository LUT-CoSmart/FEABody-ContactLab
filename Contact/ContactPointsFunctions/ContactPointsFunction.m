function [ContactPoints, ContactPointsElements] = ContactPointsFunction(ContactBody, ContactPointsName, n)
 
    nloc_cont = ContactBody.nloc;
    DofsAtNode_cont = ContactBody.DofsAtNode;
    ContactNode_cont = ContactBody.contact.nodalid;
    
    % current position of the "possible contact" nodes of the Contact Body
    ContactNodesPosition_X = ContactBody.q(xlocChosen(DofsAtNode_cont, ContactNode_cont,1)) + ...
                             ContactBody.u(xlocChosen(DofsAtNode_cont, ContactNode_cont,1));   ... % coords on X axis;
    
    ContactNodesPosition_Y =  ContactBody.q(xlocChosen(DofsAtNode_cont,ContactNode_cont,2)) + ... % coords on Y axis
                              ContactBody.u(xlocChosen(DofsAtNode_cont,ContactNode_cont,2));
    
    ContactNodesPosition = [ContactNodesPosition_X ContactNodesPosition_Y]; % nodes of the contact surfaces of the contact body     
    
    
    ContactPoints = []; % taking the first node position, otherwise later it will be omitted 
    ContactPointsElements = []; % taking the first node number

    if ~ismember(ContactPointsName, ["Gauss", "LinSpace"])
        ContactPoints = ContactNodesPosition(2:end,:); % excluding the first point
    elseif ContactPointsName == "LinSpace"
        pointfun = @geospace;
    elseif ContactPointsName == "Gauss"
        pointfun = @gauleg2;
    end  
    
    for ii = 1:size(ContactNodesPosition,1)-1
    
        % On the edge, two nodes are uniquely belong to one element only 
        ElemenNumber = find(any(nloc_cont == ii, 2) & any(nloc_cont == ii+1, 2)); 

        if ismember(ContactPointsName, ["Gauss", "LinSpace"])
            % taking two consecutive nodes
            a = ContactNodesPosition(ii,:);
            b = ContactNodesPosition(ii+1,:); 
            % Following the Intercept theorem:
            % the function devide x- and y- axes in the same propotions 
            xx = pointfun(a(1),b(1),n)'; % split in x- axis
            yy = pointfun(a(2),b(2),n)'; % split in y- axis

            % Sanity check for critical cases (needed for LinSpace)
            if abs(a(1) - b(1)) < sqrt(eps)
                xx = a(1)*ones(n,1);
            elseif abs( a(2) - b(2) ) < sqrt(eps)
                yy = a(2)*ones(n,1); 
            end

            if ContactPointsName == "LinSpace" % excluding the node from the previous devision    
                xx = xx(2:end);
                yy = yy(2:end);
            end                
            
            % To be sure of consistency
            xx = xx(:);
            yy = yy(:);

            ContactPointsElements = [ContactPointsElements; ElemenNumber*ones(length(xx),1)]; % multiplication to correlate with points
            ContactPoints = [ContactPoints; xx yy];
        else                         
            ContactPointsElements = [ContactPointsElements; ElemenNumber]; 
        end                
    end
    



