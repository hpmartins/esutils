function out = is_up(r)
    if ismember(r, [1 2 3 4 5 11 12 13 14 15])
        out = 1;
    else
        out = 0;
    end
end