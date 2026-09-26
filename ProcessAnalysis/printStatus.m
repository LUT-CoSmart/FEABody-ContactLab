function status = printStatus(approachBasis, deltaf, uu_bc, tol, ii, jj, imax, steps, titertot, Gap)
     status = false;
     if ~isnan(Gap)
        fprintf('Iteration: %d, Convergence: %10.4f, Displacements norm: %10.5f, Mean penetration: %10.7f\n', jj, norm(abs(deltaf)), norm(uu_bc), Gap);      
     else
        fprintf('Iteration: %d, Convergence: %10.4f, Displacements norm: %10.5f\n', jj, norm(abs(deltaf)), norm(uu_bc));            
     end

     if approachBasis == "Lagrange"
         cond = all(abs(uu_bc) < tol);   
     else
         cond =  all(abs(deltaf)<tol);
     end    
     
     if cond                      
         fprintf('Solution for %d / %d step  is found on %d iteration, Total CPU-time: %f\n', ii, steps, jj, titertot);         
         status = true;         
     else
         if jj==imax 
            fprintf('The solution for %d step is not found. The maximum number of iterations is reached. Total CPU-time: %.3f\n', ii, titertot);
         end   
     end 
