% function [evc, eva] = cluster_eig(H, k)
%     
%     if k > 1
%         k = nnz(diag(H)-min(diag(H)) < 30);
%     end
%     opts.disp = 1; opts.issym = 1; opts.maxit = 500;
%     fprintf('[%d]', k);
%     [evc, eva] = eigs(H, k, 'sa', opts);
%     eva = diag(eva);
% end

function [evc, eva] = cluster_eig(H, k)
    if k > 1
        m = size(H, 1);
        lwork  = 1 + 6*m + 2*(m^2);
        liwork = 3 + 5*m;
        C = lapack.lapack('SSYEVD', 'V', 'U', m, single(triu(full(H))), m, zeros(m, 1), zeros(lwork, 1), lwork, zeros(liwork, 1), liwork, 0);
        [evc, eva] = C{[4 6]};
        clear C;
    else
        opts.issym = 1;
        [evc, eva] = eigs(H, k, 'sa', opts);
        eva = diag(eva);
    end

    if k == 1
        evc = evc(:, 1:k); eva = eva(1:k);
    end
end
