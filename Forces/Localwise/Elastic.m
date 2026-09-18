function Body = Elastic(Body)

nl=Body.nElems.x * Body.nElems.y;    % computes number of elements
h = 2*sqrt(eps);
% Initialize the global stiffness
K = zeros(Body.nx,Body.nx);
Fint = zeros(Body.nx,1);

% Loop over all elements to assemble the global stiffness
for ii = 1:nl 
    % Get the vectors of coordinates and displacements for each element
    [X,U] = GetCoorDisp(ii,Body.nloc,Body.P0,Body.u);     % matrix form
    
    % Gauss integration (with separation on 2 parts) 
    % volume part 
    n_points = 2;
    [Klocv, Fev] = Stiffness(@dFe_2412V,n_points,Body.E,Body.nu,Body.Lz,U,X,h);
    % shear part  
    n_points = 1;
    [Klocs, Fes] = Stiffness(@dFe_2412S,n_points,Body.E,Body.nu,Body.Lz,U,X,h);

    % Assemble   
    Kloc=Klocs+Klocv;     
    Fintloc = Fes + Fev; % force vector of each node in element 
    % Assemble the global matrix 
    for jj = 1:8    
    	ind01 = Body.xloc(ii,jj); %Index 01
        for kk = 1:8
            ind02 = Body.xloc(ii,kk); % Index 02
            K(ind01,ind02) = K(ind01,ind02)+Kloc(jj,kk);
        end
        Fint(ind01) = Fint(ind01)+Fintloc(jj); 
    end                      
end

Body.Fint.vec = Fint;
Body.Fint.K = K;