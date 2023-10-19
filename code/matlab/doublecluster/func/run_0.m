function [b, info, h, e] = run_0(pm, wc, b_gs)

    if wc == 0
        % Base GS
        [b, info] = main_base(pm, pm.start, pm.maxholes, pm.coherent);
        
        % Hamiltoniano GS
        h = main_hamiltonian(pm, b, info);
        
        % Diagonalizacao GS
        [e.eva, e.evc] = diagonalize(h);
        e.eva = e.eva(1);
        e.evc = e.evc(:, 1);
    else
        % Base FS
        b    = cell(1, 20);
        info = cell(1, 20);
        for r = 1:20
            if (wc < 0 && pm.start(r)) || (wc > 0 && ~pm.start(r))
                new_start = pm.start;
                new_start(r) = (wc > 0);
                [b{r}, info{r}] = main_base(pm, new_start, pm.maxholes, pm.coherent);
            end
            if wc < 0 && ~pm.start(r) && pm.rs_no_elec(cc(r))
                b{r} = b_gs;
                s = other_c(r);
                b{r} = b{r}((b{r}(:, r) == 0 & b{r}(:, s) == pm.start(s)) | (pm.start(s) == 1 & b{r}(:, r) == 1 & b{r}(:, s) == 0), :);
                info{r} = calc_info(b{r}, sum(pm.start, 2), r); % AQUI: r para -Ep ou 0 para nada.
            end
        end
        
        % Hamiltonianos FS
        h = cell(1, 20);
        for r = 1:20
            if ~isempty(b{r})
                h{r} = main_hamiltonian(pm, b{r}, info{r});
                shift_noelec = one_or_two(pm.D(2)-pm.U(2), pm.D(1)-pm.U(1), cc(r)); % certo
                if wc < 0 && ~pm.start(r) && pm.rs_no_elec(cc(r)) % certo
                    h{r} = h{r} - shift_noelec*eye(size(h{r}));
                end
                if wc > 0 && pm.rs_no_elec(cc(r)) % certo
                    h{r} = h{r} + shift_noelec*eye(size(h{r}));
                end
            end
        end
        
        % Diagonalizacao FS
        e = cell(1, 20);
        for r = 1:20
            if ~isempty(h{r})
                if r > 1 && ~isempty(h{r-1}) && isequal(h{r}, h{r-1})
                    e{r}.eva = e{r-1}.eva;
                    e{r}.evc = e{r-1}.evc;
                else
                    [e{r}.eva, e{r}.evc] = diagonalize(h{r});
                end
            end
        end
    end
end
