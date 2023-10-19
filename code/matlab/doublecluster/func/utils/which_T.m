function T = which_T(pm, p)
    if ismember(p, [1 2 3 6 7 8])
        T = pm.Tp(1);
    elseif ismember(p, [4 5 9 10])
        T = pm.Ts(1);
    elseif ismember(p, [11 12 13 16 17 18])
        T = pm.Tp(2);
    elseif ismember(p, [14 15 19 20])
        T = pm.Ts(2);
    end
    
    
    if ismember(p, 1:5)
        T = pm.Tm(p)*T;
    elseif ismember(p, 6:10)
        T = pm.Tm(p-5)*T;
    elseif ismember(p, 11:15)
        T = pm.Tm(p-5)*T;
    elseif ismember(p, 16:20)
        T = pm.Tm(p-10)*T;
    end
    
end