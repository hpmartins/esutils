function out = main_hamiltonian(pm, b, info)
    ijH = zeros(100*size(b, 1), 3);
    
    nn = 0;
    
    b_lst = [sum(b > 0, 2) b];
    
    nst = min(b_lst(:, 1)) + any(sum(b == -1, 2));
    
    comb = cell(max(b_lst(:, 1)) - min(b_lst(:, 1)) + 1, 1);
    for p = min(b_lst(:, 1)):max(b_lst(:, 1))
        comb{p - min(b_lst(:, 1)) + 1} = find(b_lst(:, 1) == p);
    end

    f = 0;
    if b_lst(comb{1}(1), 1) == (nst - 1)
        for n = 1:size(comb{1}, 1)
            lst_1 = comb{1}(n);
            lst_2 = b(comb{2}, :);
            
            p = find(b(lst_1, :) == -1);
            tmp = b(lst_1, :);
            tmp(p) = 1; % Nao arrumar
            
            [~, i_d, j_d] = intersect(tmp, lst_2, 'rows');

            if size(i_d, 1) > 0
                i_d = comb{1}(i_d);
                j_d = comb{2}(j_d);
                if p <= 10
                    T = pm.Tc1;
                else
                    T = pm.Tc2;
                end

                ijH(nn+1, :) = [i_d j_d T];
                nn = nn + 1;
            end
        end
        f = 1;
    end
    
    for c = (1+f):size(comb, 1)
        if (c+1) <= size(comb, 1)
            lst_1 = b(comb{c}, :);
            lst_2 = b(comb{c+1}, :);
            
            % Ts, Tp
            for p = 1:20
                tmp_lst = lst_1;
                tmp_lst(:, p) = 1;
                [~, i_d, j_d] = intersect(tmp_lst, lst_2, 'rows');
                
                if size(i_d, 1) > 0
                    i_d = comb{c}(i_d);
                    j_d = comb{c+1}(j_d);
                    T = which_T(pm, p);
                    ijH((nn+1):(nn+size(i_d, 1)), :) = [i_d j_d T*ones(size(i_d, 1), 1)];
                    nn = nn + size(i_d, 1);
                end
            end
            
            % Tc
            for p = find(pm.coherent)
                tmp_lst = lst_1;
                tmp_lst(:, p) = 2;
                
                [~, i_d, j_d] = intersect(tmp_lst, lst_2, 'rows');
                
                if size(i_d, 1) > 0
                    i_d = comb{c}(i_d);
                    j_d = comb{c+1}(j_d);
                    
                    if p <= 10
                        T = pm.Tc1;
                    else
                        T = pm.Tc2;
                    end
                    
                    ijH((nn+1):(nn+size(i_d, 1)), :) = [i_d j_d T*ones(size(i_d, 1), 1)];
                    nn = nn + size(i_d, 1);
                end
            end
        end
    end
    
    [nEd1, nEd2, nU1, nU2, nDq1, nDq2, nJ1, nJ2, nEp, npp1, npp2, nEc1, nEc2] = calc_diag(b, info, pm);
    ii_diag = nEd1*pm.Ed1 + nEd2*pm.Ed2  + nEp*pm.Ep ...
            + nDq1*pm.Dq1 + nDq2*pm.Dq2  - nJ1*pm.J1 - nJ2*pm.J2 ...
            + nEc1*pm.Ec1 + nEc2*pm.Ec2;
        
    if pm.pp
       ii_diag = ii_diag + npp1*(pm.pps(1)-pm.ppp(1)) + npp2*(pm.pps(2)-pm.ppp(2));
    end

    if pm.kanamori
       [nu1, nul1, nu2, nul2] = calc_kanamori(b, nU1, nU2);
       ii_diag = ii_diag + nu1*pm.u1 + nul1*pm.ul1 + nu2*pm.u2 + nul2*pm.ul2;
    else
        ii_diag = ii_diag + nU1*pm.U1 + nU2*pm.U2;
    end
    
    ijH = ijH(1:(2*nn+size(b, 1)), :);
    ijH((nn+1):(2*nn), :) = [fliplr(ijH(1:nn, 1:2)) ijH(1:nn, 3)];
    ijH((2*nn+1):end,  :) = [(1:size(b, 1))' (1:size(b, 1))' ii_diag];

    out = sparse(ijH(:,1), ijH(:,2), ijH(:,3));
end