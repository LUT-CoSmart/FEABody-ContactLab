function nabla_Sigma_nn  = NablaSigma_NN(E,nu,U,X,xi,N)

        nabla_Sigma_nn = zeros(2,1);
        NablaSigma = nabla_sigma_2412(E,nu,U,X,xi(1),xi(2));
         
        i = 1;
        NablaSigma_x = [NablaSigma(1,i) NablaSigma(3,i);
                        NablaSigma(3,i) NablaSigma(2,i)];

        i = 2;
        NablaSigma_y = [NablaSigma(1,i) NablaSigma(3,i);
                        NablaSigma(3,i) NablaSigma(2,i)];

        nabla_Sigma_nn(1) = N' * NablaSigma_x * N;
        nabla_Sigma_nn(2) = N' * NablaSigma_y * N;



