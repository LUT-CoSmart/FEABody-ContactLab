function DofsFunction = AimFunctionWrapper(~, y, Body1, Body2, AimFunction)

   Body1.u=y(1:Body1.nx); 
   Body2.u=y(Body1.nx+1:end);
   DofsFunction = AimFunction(Body1,Body2);
   DofsFunction = DofsFunction(:); % numjac expects a vector
end   