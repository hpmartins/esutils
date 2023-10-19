function p = which_col(c)
    if ismember(c, [1 2 3])
        p = 1;
    elseif ismember(c, [6 7 8])
        p = 2;
    elseif ismember(c, [4 5])
        p = 3;
    elseif ismember(c, [9 10])
        p = 4;
    elseif ismember(c, [11 12 13])
        p = 5;
    elseif ismember(c, [16 17 18])
        p = 6;
    elseif ismember(c, [14 15])
        p = 7;
    elseif ismember(c, [19 20])
        p = 8;
    end
end