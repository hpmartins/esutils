N = 10;

% Valores da lista
x      = rand(N*N, 1);
index1 = reshape(repmat(1:N, N, 1), [N*N, 1]);
index2 = repmat((1:N)', N, 1);

lista = sortrows([x index1 index2], 1);

filter1 = [];
filter2 = [];
out     = [];
for i = 1:size(lista, 1)
    if ~ismember(lista(i, 2), filter1) && ~ismember(lista(i, 3), filter2)
        out(end+1, :) = lista(i, :);
        filter1(end+1) = lista(i, 2);
        filter2(end+1) = lista(i, 3);
    end
end