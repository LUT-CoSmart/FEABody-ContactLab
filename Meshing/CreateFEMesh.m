function Body = CreateFEMesh(Body)

DofsAtNode = 2;

[Body.P0,Body.nloc] = recmesh2412(Body.Lx,Body.Ly,Body.shift,Body.nElems);
Body.xloc=xlocAll2412(Body.nloc);

[nn,~] = size(Body.P0);  % get number of nodes (nn) and total number of DOFs (nx)
Body.nx = DofsAtNode*nn; % dofs (no constraints eliminated) of stucture 
u0 = zeros(Body.nx,1);   % Create the unknown vector
Body.u = u0;             % Initialize the displacements

for jj=1:nn
    q((jj-1)*DofsAtNode+1:(jj-1)*DofsAtNode+DofsAtNode)=Body.P0(jj,:); 
end  

% Define initial position
q=q(:);
Body.q = q;
   
Body.DofsAtNode = DofsAtNode;
