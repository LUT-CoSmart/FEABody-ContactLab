clc, clear, close all;

NODEs=4;
DIM=2;
POLYNOMs=4;

% Geometrical variables
syms Lx Ly Lz real;
% Material parameters
syms E nu real;
% Elements
syms xi eta real;
% nodal displacement coordinates
u = sym('u', [1 4], 'real').';
v = sym('v', [1 4], 'real').';
uu(1:2:8)=u;
uu(2:2:8)=v;
uu=uu(:);
% geometry interpolation i.e. nodal position at initial configuration
X = sym('X', [1 8], 'real').';
% Basis in reference coordinates xi={-1..1}
basis=[1,xi,eta,xi*eta];
% [Node1,Node2,Node3,Node4]
% Coordinates at nodes of element 
XI1=[-1,1,1,-1].';
XI2=[-1,-1,1,1].';
% Substituting boundary conditions
for iinode=1:NODEs
    for jjpol=1:POLYNOMs
        ii2=(iinode-1)*1+1;
        A(ii2,jjpol)=subs(basis(jjpol),[xi,eta],[XI1(iinode),XI2(iinode)]);
    end
end
% Shape functions (ansatch) in xi coordinates
Nvec=basis*A^-1;           
% Shape functions in matrix form
for ii=1:DIM
    for jj=1:POLYNOMs
        jj2=(ii-1)+(jj-1)*2+1;
        Nm(ii,jj2)=Nvec(jj);
    end
end
Nm_xi = diff(Nm, xi);
Nm_eta = diff(Nm, eta);

% Nvec_2412=matlabFunction(Nvec)
% Interpolation for assumed displacement field u 
uuh=Nm*uu;
% Interpolation for geometry (position)
XXh=Nm*X;
% Compute Je and its inverse JeInv
Je = jacobian(XXh,[xi,eta]); 
JeInv = Je^(-1);
nablau=jacobian(uuh,[xi,eta])*JeInv;
F = eye(2) + nablau;
EE=1/2*(nablau+nablau.'+nablau.'*nablau);
% strain vector
eps=[EE(1,1), EE(2,2), 2*EE(1,2)].'; 
% plain stress assumption sigmazz=0
DD=E/(1-nu^2)*[1, nu, 0;
               nu, 1, 0;
               0, 0, (1-nu)/2];
DDV=E/(1-nu^2)*[1, nu, 0;
               nu, 1, 0;
               0, 0, 0];
DDS=E/(1-nu^2)*[0, 0, 0;
               0, 0, 0;
               0, 0, (1-nu)/2];           
detJe=det(Je);% detJe is used in the FEM solver, not here.  
% Compute strain energy UdA here!                                       
UdA= 1/2*Lz*eps.'*DD*eps*detJe;        
UdAV=1/2*Lz*eps.'*DDV*eps*detJe;
UdAS=1/2*Lz*eps.'*DDS*eps*detJe;
% --------------------------------
for kk=1:8
    dFe(kk)=diff(UdA,uu(kk));
    dFeV(kk)=diff(UdAV,uu(kk));
    dFeS(kk)=diff(UdAS,uu(kk));
end 

% % matlabFunction(dFe,'file','dFe_2412','vars',{E,nu,Lz,uu,X,xi,eta});
% matlabFunction(dFeV,'file','dFe_2412V','vars',{E,nu,Lz,uu,X,xi,eta});
% matlabFunction(dFeS,'file','dFe_2412S','vars',{E,nu,Lz,uu,X,xi,eta});
% matlabFunction(Nm,'file','Nm_2412','vars',{xi,eta});
% matlabFunction(Nm_xi,'file','Nm_2412_xi','vars',{xi,eta});
% matlabFunction(Nm_eta,'file','Nm_2412_eta','vars',{xi,eta});
% matlabFunction(nablau,'file','nabla_u_2412','vars',{uu,X,xi,eta});
% 
% sigma = DD * eps;
% Sigma = [sigma(1) sigma(3);
%          sigma(3) sigma(2)];
% 
% nabla_sigma = jacobian(sigma,[xi,eta])*JeInv;
% 
% matlabFunction(nabla_sigma,'file','nabla_sigma_2412','vars',{E,nu,uu,X,xi,eta});
% matlabFunction(Sigma,'file','Sigma_2412','vars',{E,nu,uu,X,xi,eta});
% matlabFunction(F,'file','F_2412','vars',{uu,X,xi,eta});