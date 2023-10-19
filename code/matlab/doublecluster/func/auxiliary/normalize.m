function out = normalize(vec)
    vmax = max(vec);
    vmin = min(vec);
    
    if (vmax == vmin)
        out = vec;
    else
        out = (vec-vmin)/(vmax-vmin);
    end
end