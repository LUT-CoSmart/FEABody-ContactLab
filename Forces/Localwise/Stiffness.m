function [Kloc, Fe] = Stiffness(func,n_points,E,nu,Lz,U,X,h)

Fe=zeros(8,1);

[xiv,wxi]=gauleg2(-1,1,n_points);
[etav,weta]=gauleg2(-1,1,n_points);    

for ii1=1:n_points
    for jj1=1:n_points
        Fe=Fe+func(E,nu,Lz,U,X,xiv(ii1),etav(ii1))'*wxi(ii1)*weta(jj1);        
    end
end

Kloc = zeros(8,8);    % Initialize the local 
I_vec=zeros(8,1);

for jj = 1:8
    
    I_vec(jj)=1;
    Uh= U-I_vec*h;
    I_vec(jj)=0;

    Feh=zeros(8,1);
    for ii1=1:n_points
        for jj1=1:n_points
            Feh=Feh+func(E,nu,Lz,Uh,X,xiv(ii1),etav(ii1))'*wxi(ii1)*weta(jj1);
        end
    end

    Kloc(:,jj)=(Fe-Feh)/h;       
end