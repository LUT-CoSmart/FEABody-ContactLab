 
function Visualization(Body,vis, ShowNodeNumbers)
    
    FontSize = 10;
    DofsAtNode = Body.DofsAtNode;
    nl=Body.nElems.x * Body.nElems.y;    % computes number of elements
    if vis ~= "frame"
        node = reshape(Body.q + Body.u, 2, []).';
        element = Body.nloc;
        Displacementing = reshape(Body.u, 2, []).';
        texting = "Displacement: ";
        if vis == "ux"
            patching = Displacementing(:,1);
            % fprintf('presenting deformation along %s\n', vis);
        elseif vis == "uy"
            patching = Displacementing(:,2);
             % fprintf('presenting deformation along %s\n', vis);
            
        elseif vis == "u_total" 
            patching = sqrt(Displacementing(:,1).^2 + Displacementing(:,2).^2);
           %  fprintf('presenting deformation along %s\n', vis);
        elseif startsWith(vis, "sigma")
            Sigma = [];
            if vis == "sigma_VM"
                texting = "Von Mises stress: ";
                fprintf('presenting sigma_VM\n');
                
                Sigma_VM = @(S)  sqrt(S(1,1)^2  - S(1,1)*S(2,2) + S(2,2)^2 + 3*S(1,2)^2);
            else
                texting = "Inner stress: ";
                switch vis                
                    case "sigma_xx"
                        Normal_1 = [1; 0];
                        Normal_2 = Normal_1; 
                        fprintf('presenting sigma_xx\n');
                    case "sigma_yy"                
                        Normal_1 = [0; 1];
                        Normal_2 = Normal_1;
                        fprintf('presenting sigma_yy\n');
                    case "sigma_xy"
                        Normal_1 = [1; 0];
                        Normal_2 = [0; 1];
                        fprintf('presenting sigma_xy\n');
                    otherwise
                        warning('***Choose the correct axis for deformaion visualization***\n')
                        Normal_1 = [1; 0];
                        Normal_2 = Normal_1;
                        fprintf('presenting sigma_xx\n');
                end
            end
            
    
            for ii = 1:nl
                
                [X,U] = GetCoorDisp(ii,Body.nloc,Body.P0,Body.u);     % matrix form    
                S_00 = Sigma_2412(Body.E,Body.nu,U,X,-1,-1);
                S_10 = Sigma_2412(Body.E,Body.nu,U,X, 1,-1);
                S_01 = Sigma_2412(Body.E,Body.nu,U,X, 1, 1);
                S_11 = Sigma_2412(Body.E,Body.nu,U,X,-1, 1);
               
                if vis == "sigma_VM"
                    Sigma_00 = Sigma_VM(S_00);
                    Sigma_10 = Sigma_VM(S_10);
                    Sigma_01 = Sigma_VM(S_01);
                    Sigma_11 = Sigma_VM(S_11);
                else

                    Sigma_00 = Normal_1'*S_00*Normal_2;
                    Sigma_10 = Normal_1'*S_10*Normal_2;
                    Sigma_01 = Normal_1'*S_01*Normal_2;
                    Sigma_11 = Normal_1'*S_11*Normal_2;
                end
                
                Sigma = [Sigma; Sigma_00 Sigma_10 Sigma_01 Sigma_11];
            end

            
            Sigma = Sigma(:);
            [~, ia, ~] = unique(Body.nloc(:));
            patching = Sigma(ia);
            vis = strrep(vis, "sigma_", "\sigma_{") + "}"; 
        else
            warning('***Choose the correct axis for deformaion visualization***\n')
            patching = Displacementing(:,2);
            fprintf('presenting deformation along y\n', vis);
        end    
        h = patch('vertices', node, 'faces', element, 'facevertexcdata', patching, 'facecolor', 'interp');
        colormap(jet); 
        cb = colorbar;
        ylabel(cb,texting + vis, 'FontSize', [FontSize]);
        axis equal
        
        set(gca, 'FontSize', [FontSize], 'FontName','Times New Roman');
        set(text, 'FontSize', [FontSize], 'FontName','Times New Roman');
        xlabel('{\it{X}} [m]','FontName','Times New Roman','FontSize',[FontSize]),ylabel('{\it{Y}} [m]','FontName','Times New Roman','FontSize',[FontSize]),zlabel('Z [m]','FontName','Times New Roman','FontSize',[FontSize]);
    else

    end    
        
    
    
    
    
    
    
    %% Writes nodal indeces
    if ShowNodeNumbers 
        nn= Body.nx/DofsAtNode;
        [k1,k2] = size(Body.P0);
        P = zeros(k1,k2); 
        for ii=1:nn
            P(ii,:)= Body.q((ii-1)*DofsAtNode+1:(ii-1)*DofsAtNode+DofsAtNode) +...
                     Body.u((ii-1)*DofsAtNode+1:(ii-1)*DofsAtNode+DofsAtNode);    
            text(P(ii,1),P(ii,2),0, int2str(ii))
        end
    end
