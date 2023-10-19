function out = fix_L(base, info, start)
    sub = bsxfun(@minus, base, start);
    sub(sub == 2) = 0;
    
    sub = [sub sum(sub, 2)];
    
    info_col = zeros(size(base, 1), 8);
    
    for p = 1:size(sub, 1)
        if sub(p, 21) < 0
            continue;
        end
        
        cfg = sub(p, 1:20);
        
        elements = find(cfg);
        
        for el = elements
            if (el <= 10 && cfg(el) == -cfg(el+10)) || (el > 10 && cfg(el) == -cfg(el-10))
                continue;
            end
            info_col(p, which_col(el)) = info_col(p, which_col(el)) + 1;
        end
    end
    
    info(:, 12:19) = info_col;

    out = info;
end