 
function PrintResults(Body)
    
    Results = Body.results;
    disp(' Mesh & DOFs & ux & uy \\')
    for k=1:size(Results,1) 
        disp(sprintf('%dx%d & %4d  & %10.6f & %10.6f \\\\',Results(k,:)))
    end   

