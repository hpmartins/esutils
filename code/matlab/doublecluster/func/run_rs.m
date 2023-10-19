init;
%% GS
[b.gs, info.gs, h.gs, e.gs] = run_0(pm, 0, 0);
%% RS
[b.rs, info.rs, h.rs, e.rs] = run_0(pm, -1, b.gs);
%% Operator
if pm.remake
    o.rs = run_1(b.gs, b.rs, -1);
end

%% Transicões
a.rs = run_2(e.gs, o.rs, e.rs);
tmp = cat(1, a.rs{:});
shift_rs = pm.cfg.vb_p1 - min(tmp(:, 1));
for r = 1:20
    for p = 1:3
        if ~isempty(a.rs{p,r})
            a.rs{p,r}(:, 1) = a.rs{p,r}(:, 1) + shift_rs;
        end
    end
    if ~pm.start(r) && pm.rs_no_elec(cc(r))
        a.rs{2,r}(:, 2) = pm.rs_no_elec_p_wt(cc(r))*a.rs{2,r}(:, 2);
    end
end

a.all.vb = cell(2, 20);
for p = 1:2
    for r = 1:20
        if ~isempty(a.rs{p,r})
            a.all.vb{p,r} = a.rs{p,r};
        end
    end
end
a.all.vb_orb = cat(1, a.all.vb{:});
a.all.vb_orb = sortrows(a.all.vb_orb(a.all.vb_orb(:, 2) > 0.01*max(a.all.vb_orb(:, 2)), :), 1);
a.all.vb_orb = a.all.vb_orb(a.all.vb_orb(:, 1) <= 15, :);

a.vb = cell(2, 2);
for p = 1:2
    for c = 1:2
        a.vb{p,c} = cat(1, a.all.vb{p, anticc(c)});
        if isempty(a.vb{p,c})
            continue;
        end
        a.vb{p,c} = sortrows(a.vb{p,c}(a.vb{p,c}(:, 2) > 0.01*max(a.vb{p,c}(:, 2)), 1:3), 1);
        if isempty(a.vb{p,c})
            continue;
        end
        [t1, t2, ind] = consolidator(a.vb{p,c}(:, 1), a.vb{p,c}(:, 2), 'sum', 0.001);
        tmp = unique([ind a.vb{p,c}(:, 3)], 'rows');
        a.vb{p,c} = [t1 pm.weight(c)*pm.csec{p,c}*t2 tmp(:, 2)];
    end
end

%% xrange
xrange = (-5:0.1:25)';

%% Voigt
v.vb = cell(2, 2);
for p = 1:2
    if (p == 1)
        wG = pm.cfg.vb_wG_d;
    else
        wG = pm.cfg.vb_wG_p;
    end
    for c = 1:2
        v.vb{p,c} = voigtxr(a.vb{p,c}, xrange, [pm.cfg.vb_wL, pm.cfg.vb_aL, wG(c)]);
    end
end

%% Octaedro
a.oct = main_oct(pm);
a.oct(:, 2) = pm.weight(1)*pm.cfg.vb_csecp*pm.cfg.vb_oct_mult(1)*a.oct(:, 2);
a.oct(:, 3) = pm.weight(2)*pm.cfg.vb_csecp*pm.cfg.vb_oct_mult(2)*a.oct(:, 3);
v.oct1 = voigtxr(a.oct(:, [1 2]), xrange, [pm.cfg.vb_wL, pm.cfg.vb_aL, pm.cfg.vb_wG_oct(1)]);
v.oct2 = voigtxr(a.oct(:, [1 3]), xrange, [pm.cfg.vb_wL, pm.cfg.vb_aL, pm.cfg.vb_wG_oct(2)]);

%% Background
exp_data = importexp(pm.cfg.vb_exp{1});
if ~isempty(exp_data) && pm.cfg.vb_bg
    exp_data(:, 2) = normalize(exp_data(:, 2));
    exp_data = sortrows(exp_data, 1);
    exp_bg = [exp_data(:, 1) cumsum(exp_data(:, 2))];
    exp_bg(:, 2) = exp_data(end, 2)*normalize(exp_bg(:, 2));
    exp_bg = interp1(exp_bg(:, 1), exp_bg(:, 2), xrange, 'linear', 'extrap');
else
    exp_bg = zeros(size(xrange));
end

%% Normalizacao
% a.vb, v.vb, v.vb_total

max_voigts = max(v.vb{1,1} + v.vb{2,1} + v.vb{1,2} + v.vb{2,2} + v.oct1 + v.oct2);
total_a_max = max([max(a.all.vb_orb(:, 2)), max(a.oct(:, 2)), max(a.oct(:, 3))]);
for p = 1:2
    for c = 1:2
        if ~isempty(a.vb{p,c})
            v.vb{p,c} = max(v.vb{p,c})*normalize(v.vb{p,c})/max_voigts;
            a.vb{p,c}(:, 2) = 0.2*a.vb{p,c}(:, 2)/total_a_max; %nao funciona para d0 se o rs_no_elec = 0
        end
    end
end
a.oct(:, 2) = 0.2*a.oct(:, 2)/total_a_max;
v.oct1 = max(v.oct1)*normalize(v.oct1)/max_voigts;
v.oct2 = max(v.oct2)*normalize(v.oct2)/max_voigts;
v.vb_total = normalize(v.vb{1,1} + v.vb{2,1} + v.vb{1,2} + v.vb{2,2} + v.oct1 + v.oct2);
v.vb_total = normalize(v.vb_total + exp_bg);
exp_bg = v.vb_total(end)*normalize(exp_bg);


%% Plot
clf;
hold on;

% Experimental
if ~isempty(exp_data)
    scatter(exp_data(:, 1), exp_data(:, 2), '.k');
    xlim([min(exp_data(:, 1)) max(exp_data(:, 1))]);
else
    xlim([min(xrange), max(xrange)]);
end

for c = 1:2
    exp_data = importexp(pm.cfg.vb_rpes{c}); % #sorrynotsorryquevaiterquemudartodososinputs // pra plotar rpes junto, vazio pra nao plotar. no input: vb_rpes: ['exp/LMFO-RPES-Mn3d.txt', 'exp/LMFO-RPES-Fe3d.txt']
    if ~isempty(exp_data)
        exp_data(:, 2) = normalize(exp_data(:, 2));
        exp_data = sortrows(exp_data, 1);
        if (sum(exp_data(:, 1), 1) < 0)
            exp_data(:, 1) = -exp_data(:, 1);
        end
        scatter(exp_data(:, 1), 0.5*exp_data(:, 2), one_or_two('.r', '.g', c))
    end
end

if pm.cfg.vb_bg
    plot(xrange, exp_bg, ':k', 'LineSmoothing', 'on');
end
plot(xrange, v.vb_total, '-k', 'LineWidth', 2, 'LineSmoothing', 'on');
%[~, idx] = max(xrange);
%plot(xrange, v.vb_total(idx)*normalize(exp_bg), ':k', 'LineSmoothing', 'on');

plot(xrange, v.oct1, 'c', 'LineSmoothing', 'on');
plot(xrange, v.oct2, 'c', 'LineSmoothing', 'on');
barplot(a.oct(:, 1), a.oct(:,2), 2.5, 'c');
barplot(a.oct(:, 1), a.oct(:,3), 2.5, 'c');
for p = 1:2
    for c = 1:2
        if ~isempty(v.vb{p,c})
            plot(xrange, v.vb{p,c}, one_or_two(one_or_two('-r', '-g', c), one_or_two('-b', ':b', c), p), 'LineSmoothing', 'on');
            barplot(a.vb{p,c}(:, 1), a.vb{p,c}(:, 2), 2, one_or_two(one_or_two('r', 'g', c), 'b', p));
        end
    end
end
ylim([-0.1, 1.1]);
set(gca, 'XDir', 'reverse');
title(strcat(name, ' Valence Band Spectrum'));
xlabel('Binding Energy (eV)');
ylabel('Normalized Intensity');
set(gca, 'YTickLabel', '');
set(gca, 'YTick', []);
clear p r new_start t1 t2

%% Export
export(pm, sprintf('out/%s_rs_total.dat', name), [xrange, exp_bg, v.vb_total]);
export(pm, sprintf('out/%s_rs_bars_oct.dat', name), [a.oct(:, 1), a.oct(:, 2), a.oct(:, 3)]);
export(pm, sprintf('out/%s_rs_voigt_oct.dat', name), [xrange, v.oct1, v.oct2]);
for p = 1:2
    for c = 1:2
        if ~isempty(a.vb{p,c})
            export(pm, sprintf('out/%s_rs_bars_%s%d.dat', name, p2str(p), c), [a.vb{p,c}(:, 1), a.vb{p,c}(:, 2)]);
        end
        if ~isempty(v.vb{p,c})
            export(pm, sprintf('out/%s_rs_voigt_%s%d.dat', name, p2str(p), c), [xrange, v.vb{p,c}]);
        end
    end
end