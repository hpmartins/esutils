function out = infoname(info)
    s  = info(info(:, 10) == 0 & info(:, 11) == 0, :);
    
    d1s = sum(s(:, 2:4), 2);
    d2s = sum(s(:, 5:8), 2);
    
    d1 = sum(info(:, 2:4), 2);
    d2 = sum(info(:, 5:8), 2);
    
    L1 = d1 - d1s;
    L2 = d2 - d2s;
    
    C1 = info(:, 10);
    C2 = info(:, 11);
    
    out = sprintf('d%d%s%s%s%s:d%d%s%s%s%s', d1, cdt(L1 > 0, 'L', ''), cdt(L1 > 1, num2str(L1), ''), ...
                                         cdt(C1 > 0, 'C', cdt(C1 < 0, 'c', '')), cdt(C1 > 1, num2str(C1), ''), ... 
                                     d2, cdt(L2 > 0, 'L', ''), cdt(L2 > 1, num2str(L2), ''), ...
                                         cdt(C2 > 0, 'C', cdt(C2 < 0, 'c', '')), cdt(C2 > 1, num2str(C2), ''));
end

function c = cdt(condition, a , b)
    if condition
        c = a;
    else
        c = b;
    end
end