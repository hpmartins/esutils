function out = normalize2(vec)
    vmax = max(vec);
    if (vmax == min(vec))
        out = vec;
    else
        out = vec/vmax;
    end
end