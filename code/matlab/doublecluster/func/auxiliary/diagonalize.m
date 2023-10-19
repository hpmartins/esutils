function [eva, evc] = diagonalize(h)
    if (size(h, 1) == 1)
        evc = 1;
        eva = full(h);
    else
        [evc, eva] = cluster_eig(h, size(h, 1));
        fprintf(sprintf('%s\n', any(evc == 1.0)));
    end
end