function Body = SetMaterial(Body,name)
    
    if nargin < 2
        warning('unknown Material "%s". Use "KS" or "Neo", set to "KS"', name);
        name = "KS";
    end

    name = string(name);               
    if ~ismember(name, ["KS","Neo"])
        warning('unknown Material "%s". Use "KS" or "Neo", set to "KS"', name);
        name = "KS";
    end

    Body.ForceVolumeFunction = str2func("dFe_2412V_" + name);        
    Body.ForceShearFunction  = str2func("dFe_2412S_" + name);        
    Body.StressFunction      = str2func("Sigma_raw_2412_" + name); 