function [op_d, op_p, op_c] = main_operator(is, fs, orb, signal)
    % d
    if signal > 0
        old = 0;
        new = 1;
    else
        old = 1;
        new = 0;
    end
    is_lst = find(is(:, orb) == old);
    is_alt = is(is(:, orb) == old, :);
    is_alt(:, orb) = new;
    [~, i_d, j_d] = intersect(is_alt, fs, 'rows');
    i_d = is_lst(i_d);
    
    % p
    [~, i_p, j_p] = intersect(is, fs, 'rows');
    
    % c (vem do banho)
    if signal > 0
        old = 0;
        new = 2;
    else
        old = 2;
        new = 0;
    end
    is_lst = find(is(:, orb) == old);
    is_alt = is(is(:, orb) == old, :);
    is_alt(:, orb) = new;
    [~, i_c, j_c] = intersect(is_alt, fs, 'rows');
    i_c = is_lst(i_c);

    % c (vai pro banho)
    if signal < 0
        old = -1;
        new = 0;
        is_lst = find(is(:, orb) == old);
        is_alt = is(is(:, orb) == old, :);
        is_alt(:, orb) = new;
        [~, i_c2, j_c2] = intersect(is_alt, fs, 'rows');
        i_c2 = is_lst(i_c2);
        
        i_c = [i_c; i_c2];
        j_c = [j_c; j_c2];
    end
    
    nn_d = size(i_d, 1);
    nn_p = size(i_p, 1);
    nn_c = size(i_c, 1);
    
    op_d = sparse(i_d, j_d, ones(nn_d, 1), size(is, 1), size(fs, 1));
    op_p = sparse(i_p, j_p, ones(nn_p, 1), size(is, 1), size(fs, 1));
    op_c = sparse(i_c, j_c, ones(nn_c, 1), size(is, 1), size(fs, 1));
end