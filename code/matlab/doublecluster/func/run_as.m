init;
%% GS
[b.gs, info.gs, h.gs, e.gs] = run_0(pm, 0, 0);
%% AS
[b.as, info.as, h.as, e.as] = run_0(pm, 1, 0);
%% Operator
if pm.remake
    o.as = run_1(b.gs, b.as, 1);
end

%% Transicões
a.as = run_2(e.gs, o.as, e.as);
a.all.as = cell(3, 2);
for p = 1:3
    for c = 1:2
        a.all.as{p,c} = cat(1, a.as{p, anticc(c)});
        if isempty(a.all.as{p,c})
            continue;
        end
        a.all.as{p,c} = sortrows(a.all.as{p,c}(a.all.as{p,c}(:, 2) > 0.0, 1:3), 1);
        if isempty(a.all.as{p,c})
            continue
        end
        [t1, t2, ind] = consolidator(a.all.as{p,c}(:, 1), a.all.as{p,c}(:, 2), 'sum', 0.01);
        tmp = unique([ind a.all.as{p,c}(:, 3)], 'rows');
        a.all.as{p,c} = [t1 pm.weight(c)*t2 tmp(:, 2)]; %todos os inputs de STRO tem essa flag 
    end                                                    %logo no início, caso queira copiar..
end
a.all.as_orb = cat(1, a.as{:});
a.all.as_orb = sortrows(a.all.as_orb(a.all.as_orb(:, 2) > 0.001*max(a.all.as_orb(:, 2)), :), 1);
a.all.as_orb = a.all.as_orb(a.all.as_orb(:, 1) <= 15, :);

for c = 1:2
    a.all.as{2,c} = a.all.as{2,c}(a.all.as{2,c}(:,2)/max([a.all.as{2,1}(:, 2); a.all.as{2,2}(:, 2)]) >= 1e-3, :);
end
%% Voigt
tmp = cat(1, a.all.as{2, :});
shift_xas = pm.cfg.xas_p1 - min(tmp(:, 1));
for c = 1:2
    a.all.as{2,c}(:, 1) = a.all.as{2,c}(:, 1) + shift_xas;
end
xrange = transpose(linspace(min(tmp(:,1)) - 10, max(tmp(:,1)) + 10, 500)) + shift_xas;

v.as = cell(1, 2);
for c = 1:2
    v.as{c} = zeros(size(xrange));
    if ~isempty(a.all.as{2,c})
        v.as{c} = voigtxr(a.all.as{2,c}(:, 1:2), xrange, [pm.cfg.xas_wL, pm.cfg.xas_aL, pm.cfg.xas_wG]);
    end
end
v.all.as = v.as{1} + v.as{2};
v.as{1} = v.as{1}/max(v.all.as);
v.as{2} = v.as{2}/max(v.all.as);
v.all.as = v.all.as/max(v.all.as);

%% Export
export(pm, sprintf('out/%s_xas_total.dat', name), [xrange, v.as{1}, v.as{2}, v.all.as]);
for c = 1:2
    if ~isempty(a.all.as{2,c})
        export(pm, sprintf('out/%s_xas_bars_p%d.dat', name, c), [a.all.as{2,c}(:, 1), a.all.as{2,c}(:, 2)]);
    end
end

%% Plot
clf;
hold on;

exp_data = importexp(pm.cfg.xas_exp{1});
if ~isempty(exp_data)
    scatter(exp_data(:, 1), exp_data(:, 2), '.k');
end
if (strcmp(name, 'SFMO'))
    xas_a_mult = 5;
    xas_v_mult = 2*0.70/max(v.all.as);
else
    xas_a_mult = 1;
    xas_v_mult = 1;
end
plot(xrange, xas_v_mult*v.all.as/2, '-k', 'LineSmoothing', 'on', 'LineWidth', 1.5);
for c = 1:2
    if ~isempty(a.all.as{2,c})
        plot(xrange, v.as{c}/2, one_or_two('r', 'g', c), 'LineSmoothing', 'on');
        barplot(a.all.as{2,c}(:, 1), xas_a_mult*a.all.as{2,c}(:, 2), 2, one_or_two('r', 'g', c));
    end
end
%xlim([min(exp_data(:, 1)) max(exp_data(:, 1))]);
xlim([min(exp_data(:, 1)) max(exp_data(:, 1))]);
if (strcmp(name, 'SFMO'))
    xlim([525 538]);
end
ylim([-0.1, 1.1]);
title(strcat(name, ' O 1s X-ray Absorption Spectrum'));
xlabel('Energy (eV)');
ylabel('Normalized Intensity');
set(gca, 'YTickLabel', '');
set(gca, 'YTick', []);


clear p r new_start t1 t2