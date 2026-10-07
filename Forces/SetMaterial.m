function Body = SetMaterial(Body,name)
    

    
    if (nargin < 2) || (name ~= "KS" || name ~= "Neo")
        warning('unknown Material "%s". Use "KS" or "Neo", set to "KS"', name);
        name = "KS";
    end

    path = 'ElementFunctions\' + name;
    PackageFolder = fullfile(path);
    fullfile(PackageFolder, @dFe_2412V); 
     
    Body.ForceVolumeFunction = fullfile(PackageFolder, @dFe_2412V) ;               % normal-strain part, 2x2 Gauss points
     Body.ForceShearFunction  =  fullfile(PackageFolder, @dFe_2412S);               % shear part, 1 Gauss point
     Body.StressFunction      = fullfile(PackageFolder, @Sigma_raw_2412);         % 2nd Piola-Kirchhoff stress
    