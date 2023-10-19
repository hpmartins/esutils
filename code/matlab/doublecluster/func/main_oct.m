function oct = main_oct(pm)
    pos = [3 1  1   0  0  0;
           1 1  2  -2  0  0;
           1 1  3  -1  1  0;
           1 1  4  -1 -1  0;
           1 1  5  -1  0  1;
           1 1  6  -1  0 -1;
           2 1  7   1 -1  0;
           2 1  8   1  1  0;
           2 1  9   2  0  0;
           2 1 10   1  0 -1;
           2 1 11   1  0  1;
           1 2 12  -1  0  0;
           2 2 13   1  0  0];
%         pos = [1 1  1   0  0  0;
%                1 1  2  -2  0  0;
%                2 1  2   2  0  0;
%                1 1  3  -1  1  0;
%                2 1  3   1 -1  0;
%                1 1  4  -1 -1  0;
%                2 1  4   1  1  0;
%                1 1  5  -1  0  1;
%                2 1  5   1  0 -1;
%                1 1  6  -1  0 -1;
%                2 1  6   1  0  1;
%                1 2  7  -1  0  0;
%                2 2  8   1  0  0];
    
    tmp = unique(pos(:, [2 3]), 'rows');
    H_size = 3*sum(tmp(:,1) == 1) + 5*sum(tmp(:,1) == 2);
    H_oct = zeros(H_size);

    for p = 1:size(pos, 1)
        for q = p:size(pos, 1)
            t1 = pos(p,2); t2 = pos(q,2);
            i  = pos(p,3); j  = pos(q,3);
            
            c = pos(q, 1);
            % diagonal
            if i == j
                if t1 == 1
                    for m = 0:2
                        H_oct(nst(i)+m, nst(i)+m) = pm.Ep;
                    end
                end
                if t1 == 2
                    for m = 0:4
                        Ed = pm.Ed(c) + sum(pm.start(1, one_or_two(1:10, 11:20, c)),2)*pm.U(c);
                        H_oct(nst(i)+m, nst(i)+m) = Ed + 6*pm.Dq(c)*delta(m <= 1, 1) - 4*pm.Dq(c)*delta(m >= 2, 1);
                    end
                end
            else
                r = pos(q,end-2:end) - pos(p,end-2:end);
                if norm(r) < 1.5 && norm(r) > 0
                    for m = 0:one_or_two(2, 4, t1)
                        for n = 0:one_or_two(2, 4, t2)
                            T = slater(pm, r, c, t1, t2, m+1, n+1);
                            pos1 = nst(i) + m;
                            pos2 = nst(j) + n;
                            H_oct(pos1, pos2) = T;
                            H_oct(pos2, pos1) = T;
                        end
                    end
                end
            end
        end
    end

    [evc_oct, eva_oct] = eig(H_oct); eva_oct = diag(eva_oct);
    char_p1 = zeros(H_size, 1);
    char_p2 = zeros(H_size, 1);
    char_d1 = zeros(H_size, 1);
    char_d2 = zeros(H_size, 1);

    d1_l = d_lines(1);
    d2_l = d_lines(2);
    p1_l = p_lines(1);
    p2_l = p_lines(2);
    pc_l = p_lines(3);

     for q = 1:H_size
         vec = evc_oct(:, q);
         char_p1(q) = sum(vec(p1_l).^2, 1) + sum(vec(pc_l).^2, 1)/2;
         char_p2(q) = sum(vec(p2_l).^2, 1) + sum(vec(pc_l).^2, 1)/2;
         char_d1(q) = sum(vec(d1_l).^2, 1);
         char_d2(q) = sum(vec(d2_l).^2, 1);
     end
     
     tmp2 = [-eva_oct 2*char_p1 2*char_p2 2*char_d1 2*char_d2];
     [E, p1, ind] = consolidator(roundn(tmp2(:, 1), -3), roundn(tmp2(:, 2), -3), 'sum', 0.01);
     [~, p2] = consolidator(ind, roundn(tmp2(:, 3), -3), 'sum', 0.01);
     [~, d1] = consolidator(ind, roundn(tmp2(:, 4), -3), 'sum', 0.01);
     [~, d2] = consolidator(ind, roundn(tmp2(:, 5), -3), 'sum', 0.01);
     oct = sortrows([E, p1, d1, p2, d2, (p1+d1+p2+d2), (d1+d2)./(p1+p2)], 1);
     tmp = sortrows(oct, 7);
     oct = sortrows(tmp(cumsum(tmp(:, 6)) <= 46.0, [1 2 4 3 5]), 1);
     
    function out = p_lines(cc)
        list = pos((pos(:,2) == 1) & (pos(:,1) == cc), 3);
        out = [];
        for kk = 1:numel(list)
            out = [out, nst(list(kk)):(nst(list(kk))+2)];
        end
        out = unique(out);
    end
    function out = d_lines(cc)
        list = pos((pos(:,2) == 2) & (pos(:,1) == cc), 3);
        out = [];
        for kk = 1:numel(list)
            out = [out, nst(list(kk)):(nst(list(kk))+4)];
        end
        out = unique(out);
    end
    function out = nst(nn)
        out = 1;
        teste = unique(pos(:, [2 3]), 'rows');
        for kk = 1:(nn-1)
            out = out + one_or_two(3, 5, teste(kk,1));
            if kk == nn
                break;
            end
        end
    end
end