function out = calc_info(b, nst, r)
    if cc(r) == 1
        is_c1 = 1;
    else
        is_c1 = 0;
    end
    d1 = sum(b(:, 1:10) > 0, 2);
    d2 = sum(b(:, 11:20) > 0, 2);
    L  = sum(b == 1, 2) - nst + (r > 0) + sum(b == -1, 2); %tirei o 2 pra simular 'direito' outras coisas
% +   (r > 0): d0L   - Ti no topo
%   0*(r > 0): d0    - Ti no topo, resto shifta pra esquerda
% -   (r > 0): d0-L  - Ti no topo, resto shifta o dobro pra esquerda
% + 2*(r > 0): d0L2  - Perfeito
    t2g_1 = sum(b(:, [1 2 3 6 7 8]) == 1, 2) + (r > 0)*is_c1*is_t2g(r); % nao conta coerentes
    eg_1 = sum(b(:, [4 5 9 10]) == 1, 2) + (r > 0)*is_c1*(~is_t2g(r));
    t2g_2 = sum(b(:, [11 12 13 16 17 18]) == 1, 2) + (r > 0)*(~is_c1)*is_t2g(r);
    eg_2 = sum(b(:, [14 15 19 20]) == 1, 2) + (r > 0)*(~is_c1)*(~is_t2g(r));
    out = [d1 d2 L t2g_1 eg_1 t2g_2 eg_2 r*ones(size(d1, 1), 1)];
end