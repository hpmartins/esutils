function [out, info] = main_base(pm, start, maxholes, coherent)
    all_cfgs = npermutek([0 1], 20);
    nst = sum(start, 2);

    % Base com os mesmos eletrons do start
    fs = find(start);
    base_def = all_cfgs(sum(all_cfgs(:, fs), 2) == numel(fs) & sum(all_cfgs, 2) <= nst+maxholes, :);
 
    % Base coerente
    base_coherent = [];
    if any(coherent)
        base_coherent = zeros(size(base_def, 1), 20);
        nn = 0;
        for i = find(coherent)
            if start(i) == 1
                coherent_state = start;
                coherent_state(i) = -1;
                base_coherent(nn+1, :) = coherent_state;
                nn = nn + 1;
            end
            coherent_list = base_def(base_def(:, i) == 1 & start(i) == 0, :);
            coherent_list(:, i) = 2;
            base_coherent((nn+1):(nn+size(coherent_list, 1)), :) = coherent_list;
            nn = nn + size(coherent_list, 1);
        end
        base_coherent = base_coherent(1:nn, :);
    end
    
    % Base com transferencia de um cluster pro outro
    if pm.bs_transfer
        base_transfer = zeros(20*size(base_def, 1), 20);
        nn = 0;
        for i = find(start)
            j = other_c(i);
            possible_transfers = base_def(base_def(:, j) == 0 & sum(base_def, 2) <= nst, :);
            possible_transfers(:, i) = 0;
            possible_transfers(:, j) = 1;
            base_transfer((nn+1):(nn+size(possible_transfers, 1)), :) = possible_transfers;
            nn = nn + size(possible_transfers, 1);

            if any(coherent)
                possible_transfers_coh = base_coherent(base_coherent(:, j) == 0, :);
                possible_transfers_coh(any(possible_transfers_coh == -1, 2), :) = [];
                possible_transfers_coh(:, i) = 0;
                possible_transfers_coh(:, j) = 1;
                base_transfer((nn+1):(nn+size(possible_transfers_coh, 1)), :) = possible_transfers_coh;
                nn = nn + size(possible_transfers_coh, 1);
            end
        end
        base_transfer = base_transfer(1:nn, :);
    else
        base_transfer = [];
    end

    out = unique([base_def; base_coherent; base_transfer], 'rows', 'stable');
    
    tmp = [sum(out, 2) out];
    out = sortrows(tmp, 1);
    out = out(:, 2:end);
    
    info = calc_info(out, nst, 0);
end
