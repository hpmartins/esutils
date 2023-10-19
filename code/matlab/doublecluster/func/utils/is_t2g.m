function out = is_t2g(r)
    if ismember(r, [1 2 3 6 7 8 11 12 13 16 17 18])
        out = 1;
    else
        out = 0;
    end
end