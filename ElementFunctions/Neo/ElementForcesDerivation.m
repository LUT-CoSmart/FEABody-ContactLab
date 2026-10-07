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
detJe=det(Je);
mu = E/(2*(1+nu));
lambda   = E*nu/((1+nu)*(1-2*nu));
C = F.'*F;             
C_inv = inv(C);
J = det(F);                             
W = mu/2*(trace(C) + 1 - 3) - mu*log(J) + lambda/2*log(J)^2;
UdAV = Lz*W*detJe; 
UdAS = sym(0);          % no separate shear part -> dFe_2412S returns zeros
Sigma = mu*(eye(2) - C_inv) + lambda*log(J)*C_inv;

% --------------------------------
for kk=1:8
    dFeV(kk)=diff(UdAV,uu(kk));
    dFeS(kk)=diff(UdAS,uu(kk));
end

matlabFunction(dFeV, 'file','dFe_2412V',     'vars',{E,nu,Lz,uu,X,xi,eta});
matlabFunction(dFeS, 'file','dFe_2412S',     'vars',{E,nu,Lz,uu,X,xi,eta});
matlabFunction(Sigma,'file','Sigma_raw_2412','vars',{E,nu,uu,X,xi,eta});