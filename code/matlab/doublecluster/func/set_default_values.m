function npm = set_default_values(opm)
    npm = opm;
    
    if ~(isfield(npm, 'coherent'))
        npm.coherent = [0,0,0,0,0, 0,0,0,0,0,   0,0,0,0,0, 0,0,0,0,0];
    end
    
    if ~isfield(npm, 'rs_no_elec')
        npm.rs_no_elec = [0, 0];
    end
    
    if ~isfield(npm, 'rs_no_elec_p_wt')
        npm.rs_no_elec_p_wt = [0, 0];
    end
    
    if ~isfield(npm, 'weight')
        npm.weight = [1, 1];
    end
    
    if ~isfield(npm, 'bs_transfer')
        npm.bs_transfer = 0;
    end
    
    if ~isfield(npm, 'Tm')
        npm.Tm = [1, 1, 1, 1, 1,     1, 1, 1, 1, 1];
    end
    
    if ~isfield(npm, 'Dc1')
        npm.Dc1 = 0;
    end
    if ~isfield(npm, 'Tc1')
        npm.Tc1 = 0;
    end
    if ~isfield(npm, 'Dc2')
        npm.Dc2 = 0;
    end
    if ~isfield(npm, 'Tc2')
        npm.Tc2 = 0;
    end
end