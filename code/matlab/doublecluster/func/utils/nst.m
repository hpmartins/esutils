function out = nst(pos, nn)
    out = 1;
    teste = unique(pos(:, [2 3]), 'rows');
    for kk = 1:(nn-1)
        out = out + one_or_two(3, 5, teste(kk,1));
        if kk == nn
            break;
        end
    end
end