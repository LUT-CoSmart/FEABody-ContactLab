function PostProcess(Body1, Body2, ShowVisualization, WhatoToShow, ShowNodeNumbers, approach, ContactPointfunc,Gapfunc,Solution,PointsofInterest)

    if approachSubtype == "penalty-Nitsche"
        alpha_eff = alpha * approach.lambda.meaning / approach.penalty;
        fprintf('penalty-Nitsche: alpha = %g, alpha_eff: mean = %g, max = %g\n',alpha, mean(alpha_eff), max(alpha_eff));
    end
    
    PrintResults(Body1)
    PrintResults(Body2)
    
    if ShowVisualization
        hold on 
        Visualization(Body1,WhatoToShow,ShowNodeNumbers);
        Visualization(Body2,WhatoToShow,ShowNodeNumbers);
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % visualization of contact points and contact method
        [ContactPoints, ~] = ContactPointfunc(Body1);  
        h1 = plot(ContactPoints(:,1),ContactPoints(:,2),'ok','MarkerFaceColor', 'k', 'MarkerSize', 10);
        legend('contact points')
        legend(h1, 'contact points');
        Gap =  Gapfunc(Body1,Body2);
        gapStr = sprintf('%.5f', Gap);
        fullstr = "Solution method - " + Solution +", Contact method = " + approach.Name + ...
                   ", Contact points choice = " + PointsofInterest.Name + ", Total Gap = " + gapStr;
        title(fullstr, 'Interpreter', 'latex','FontSize', 15);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end    
