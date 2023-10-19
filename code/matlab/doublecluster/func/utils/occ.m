function occ(evc, dc)
    if size(evc, 2) > 1
        fprintf('Erro! Primeiro argumento (autovetor) deve ser apenas uma coluna.\n');
        return;
    end
    
    pesos = evc.*evc;
    occup_1   = sum(pesos.*sum(dc(:, 1:10), 2));
    occup_2   = sum(pesos.*sum(dc(:, 11:end), 2));
    

    lista = [info_base(dc) 100*pesos];
    lista(lista(:, 8) < 1e-10, 8) = 0;
    

    
    [un_lista, ~, un_ind] = unique(lista(:, 1:7), 'rows');
    lista = sortrows([un_lista accumarray(un_ind, lista(:, 8))], -8);
    
        
    
    to_show = lista(lista(:, 8) > -0.1, :);
    
    for i = 1:size(to_show, 1)
        fprintf('[%10s] %5.2f\n', cfgname(to_show(i, 1:7)), to_show(i, 8));
    end
    fprintf(' <nd> = %.3f / %.3f\n', occup_1, occup_2);
    
    % Ocupa��o por cluster
%     for c = 1:2
%         if c == 1
%             range = 1:3;
%         else
%             range = 4:6;
%         end
%         [un_lista, ~, un_ind] = unique(to_show(:, range), 'rows');
%         tmp = sortrows([un_lista accumarray(un_ind, to_show(:, 7))], -4);
%         
%         fprintf(' C%d -- ', c);
%         for i = 1:size(tmp, 1)
%             fprintf('%s: %.2f%%', cfgname(tmp(i, [1 2 3])), tmp(i, 4));
%             if i ~= size(tmp, 1)
%                 fprintf(', ');
%             end
%         end
%         fprintf('\n');
%     end
end


function c = cdt(condition, a , b)
    if condition
        c = a;
    else
        c = b;
    end
end

function str = cfgname(q)
    str = sprintf('d%d%s%s%s%s%s%s:d%d%s%s%s%s', q(1), ...
                             cdt(q(2) > 0, 'C', ''), cdt(q(2) > 1, num2str(q(2)), ''), ...
                             cdt(q(3) > 0, 'c', ''), cdt(q(3) > 1, num2str(q(3)), ''), ...
                             cdt(q(4) > 0, ':L', ''), cdt(q(4) > 1, num2str(q(4)), ''), ...
                             q(5), ...
                             cdt(q(6) > 0, 'C', ''), cdt(q(6) > 1, num2str(q(6)), ''), ...
                             cdt(q(7) > 0, 'c', ''), cdt(q(7) > 1, num2str(q(7)), ''));
end