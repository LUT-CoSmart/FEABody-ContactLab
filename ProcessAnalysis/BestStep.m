function [Body1,Body2,Gap] = BestStep(bestStep,lambda,Gap,Body1,Body2,deltaf,ii,imax)
       
    persistent deltaf_best Body1_best Body2_best Gap_best iteration

    if bestStep.meaning   
               
       delta =  norm(abs(deltaf)); 
       
       if lambda == 1 && ii == 1 
       % lambda here is accounted only if bestStep activated and we are already on another step of backtracking
       % thus searching and taking the best step not within current lambda but all; 
       
          iteration = ii;
          deltaf_best = delta;
          Body1_best = Body1;
          Body2_best = Body2;
          Gap_best = Gap; 

       else 

          if delta < deltaf_best
             deltaf_best = delta;
             Body1_best = Body1;
             Body2_best = Body2;
             Gap_best = Gap;
             iteration = ii;
          end

       end
    
       if ii==imax               
          Body1 = Body1_best;
          Body2 = Body2_best;
          Gap =  Gap_best;    
          warning('Chosen the following solution: Iteration: %d, Convergence: %10.4f, Total gap: %10.7f\n', iteration, deltaf_best, Gap); 
       end 

    end  