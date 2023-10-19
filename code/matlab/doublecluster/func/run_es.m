init;
%% GS
[b.gs, info.gs, h.gs, e.gs] = run_0(pm, 0, 0);
%% AS
[b.as, info.as, h.as, e.as] = run_0(pm, 1, 0);
%% RS
[b.rs, info.rs, h.rs, e.rs] = run_0(pm, -1, b.gs);
%% Operator
if pm.remake
    o.as = run_1(b.gs, b.as, 1);
    o.rs = run_1(b.gs, b.rs, -1);
end

%% Transicões AS
a.as = run_2(e.gs, o.as, e.as);
a.all.as = cell(3, 2);
for p = 1:3
    for c = 1:2
        a.all.as{p,c} = cat(1, a.as{p, anticc(c)});
        if isempty(a.all.as{p,c})
            continue;
        end
        a.all.as{p,c} = sortrows(a.all.as{p,c}(a.all.as{p,c}(:, 2) > 0.0001*max(a.all.as{p,c}(:, 2)), 1:3), 1);
        if isempty(a.all.as{p,c})
            continue
        end
        [t1, t2, ind] = consolidator(a.all.as{p,c}(:, 1), a.all.as{p,c}(:, 2), 'sum', 0.001);
        tmp = unique([ind a.all.as{p,c}(:, 3)], 'rows');
        a.all.as{p,c} = [t1 pm.weight(c)*t2 tmp(:, 2)];
    end
end
a.all.as_orb = cat(1, a.as{:});
a.all.as_orb = sortrows(a.all.as_orb(a.all.as_orb(:, 2) > 0.001*max(a.all.as_orb(:, 2)), :), 1);
a.all.as_orb = a.all.as_orb(a.all.as_orb(:, 1) <= 15, :);

%% Transicões RS
a.rs = run_2(e.gs, o.rs, e.rs);
for r = 1:20
    if ~pm.start(r) && pm.rs_no_elec(cc(r))
        a.rs{2,r}(:, 2) = pm.rs_no_elec_p_wt(cc(r))*a.rs{2,r}(:, 2);
    end
end
a.all.rs = cell(3, 2);
for p = 1:3
    for c = 1:2
        a.all.rs{p,c} = cat(1, a.rs{p, anticc(c)});
        if isempty(a.all.rs{p,c})
            continue;
        end
        a.all.rs{p,c} = sortrows(a.all.rs{p,c}(a.all.rs{p,c}(:, 2) > 0.001*max(a.all.rs{p,c}(:, 2)), 1:3), 1);
        if isempty(a.all.rs{p,c})
            continue;
        end
        [t1, t2, ind] = consolidator(a.all.rs{p,c}(:, 1), a.all.rs{p,c}(:, 2), 'sum', 0.001);
        tmp = unique([ind a.all.rs{p,c}(:, 3)], 'rows');
        a.all.rs{p,c} = [t1 pm.weight(c)*t2 tmp(:, 2)];
    end
end
a.all.rs_orb = cat(1, a.rs{:});
a.all.rs_orb = sortrows(a.all.rs_orb(a.all.rs_orb(:, 2) > 0.001*max(a.all.rs_orb(:, 2)), :), 1);
a.all.rs_orb = a.all.rs_orb(a.all.rs_orb(:, 1) <= 15, :);

shift_es = min(a.all.rs_orb(:, 1)) - pm.cfg.vb_p1;

a.oct = main_oct(pm);
a.oct(:, 1) = -a.oct(:, 1);
a.oct(:, 2) = pm.weight(1)*pm.cfg.vb_oct_mult(1)*a.oct(:, 2);
a.oct(:, 3) = pm.weight(2)*pm.cfg.vb_oct_mult(2)*a.oct(:, 3);

%% Junção AS-RS
a.es = cell(2, 4);
for r = 1:20
    rs_d = cat(1, a.rs{[1,3], r});
    rs_p = a.rs{2,r};
    as_d = cat(1, a.as{[1,3], r});
    as_p = a.as{2,r};
    
    col  = (1 + ~is_up(r)) + 2*(cc(r) - 1);
    spin = one_or_two(-1, 1, 1+is_up(r));
    
    if ~isempty(rs_d)
        rs_d = [-rs_d(:, 1), rs_d(:, 2)];
    end
    
    if ~isempty(as_d)
        as_d = as_d(:, 1:2);
    end
    
    if ~isempty(rs_p)
        rs_p = [-rs_p(:, 1), rs_p(:, 2)];
    end
    
    if ~isempty(as_p)
        as_p = as_p(:, 1:2);
    end

    a.es{1,col} = [a.es{1,col}; rs_d -1*ones(size(rs_d, 1), 1); as_d ones(size(as_d, 1), 1)];
    a.es{2,col} = [a.es{2,col}; rs_p -1*ones(size(rs_p, 1), 1); as_p ones(size(as_p, 1), 1)];
end

for i = 1:2
    for j = 1:4
        if ~isempty(a.es{i,j})
            [t1, t2, ind] = consolidator(a.es{i,j}(:, 1), a.es{i,j}(:, 2), 'sum', 0.001);
            t3 = unique([ind a.es{i,j}(:, 3)], 'rows');
            a.es{i,j} = [t1+shift_es t2 t3(:, 2)];
        end
    end
end

%% Voigt
v.es = cell(2, 4);
for i = 1:2
    for j = 1:4
        if ~isempty(a.es{i,j})
            v.es{i,j} = voigt(a.es{i,j}(:, 1:2), [pm.cfg.vb_wL, pm.cfg.vb_aL, pm.cfg.vb_wG_d(c)], [-15 15 5 0.05]);
        end
    end
end

v.oct1 = voigtxy(a.oct(:, 1), a.oct(:, 2)/2, [pm.cfg.vb_wL, pm.cfg.vb_aL, pm.cfg.vb_wG_oct(c)], [-15 15 5 0.05]);
v.oct2 = voigtxy(a.oct(:, 1), a.oct(:, 3)/2, [pm.cfg.vb_wL, pm.cfg.vb_aL, pm.cfg.vb_wG_oct(c)], [-15 15 5 0.05]);

v.es_total = cell(1, 2);
for i = 1:2
    tmp = cat(1, v.es{:, [i, i+2]});
    [t1, t2, ~] = consolidator(tmp(:, 1), tmp(:, 2), 'sum', 0.001);
    v.es_total{i} = [t1 t2];
    v.es_total{i}(:, 2) = v.es_total{i}(:, 2) + v.oct1(:, 2) + v.oct2(:, 2);
end

%% Export

%% Plot
clf;
hold on;

clr.d_rs1 = [1 0 0];
clr.d_rs2 = [0 1 0];
clr.p_rs1 = [0 0 1];
clr.p_rs2 = [0 0 1];

clr.d_as1 = [0.8 0.1 0.0];
clr.d_as2 = [0.0 0.5 0.0];
clr.p_as1 = [0.0 0.5 0.8];
clr.p_as2 = [0.0 0.5 0.8];
clr.oct   = [0.0 0.8 1.0];

clrs = {{clr.d_rs1, clr.d_as1; clr.p_rs1, clr.p_as1}, {clr.d_rs2, clr.d_as2; clr.p_rs2, clr.p_as2}};

for i = 1:2
    for j = 1:4
        if ~isempty(a.es{i,j})
            rs_b = a.es{i,j}(a.es{i,j}(:, 3) < 0, :);
            as_b = a.es{i,j}(a.es{i,j}(:, 3) > 0, :);
            
            rs_c = one_or_two(clrs{1}(:, 1), clrs{2}(:, 1), 1 + (j >= 3));
            as_c = one_or_two(clrs{1}(:, 2), clrs{2}(:, 2), 1 + (j >= 3));
            
            barplot(rs_b(:, 1), one_or_two(-1, 1, 1+mod(j,2))*rs_b(:, 2)/3, 2, rs_c{i});
            barplot(as_b(:, 1), one_or_two(-1, 1, 1+mod(j,2))*as_b(:, 2)/3, 2, as_c{i});
            
            plot(v.es{i,j}(:, 1), one_or_two(-1, 1, 1+mod(j,2))*v.es{i,j}(:, 2), 'Color', rs_c{i}, 'LineSmoothing', 'on');
            
            if (j <= 2)
                export(pm, sprintf('out/%s_es_bars_%s1_%s.dat', name, one_or_two('d', 'p', i), one_or_two('dn', 'up', 1+mod(j,2))), [a.es{i,j}(:, 1), one_or_two(-1, 1, 1+mod(j,2))*a.es{i,j}(:,2)]);
                export(pm, sprintf('out/%s_es_voigt_%s1_%s.dat', name, one_or_two('d', 'p', i), one_or_two('dn', 'up', 1+mod(j,2))), [v.es{i,j}(:, 1), one_or_two(-1, 1, 1+mod(j,2))*v.es{i,j}(:,2)]);
            else
                export(pm, sprintf('out/%s_es_bars_%s2_%s.dat', name, one_or_two('d', 'p', i), one_or_two('dn', 'up', 1+mod(j,2))), [a.es{i,j}(:, 1), one_or_two(-1, 1, 1+mod(j,2))*a.es{i,j}(:,2)]);
                export(pm, sprintf('out/%s_es_voigt_%s2_%s.dat', name, one_or_two('d', 'p', i), one_or_two('dn', 'up', 1+mod(j,2))), [v.es{i,j}(:, 1), one_or_two(-1, 1, 1+mod(j,2))*v.es{i,j}(:,2)]);
            end
        end
    end
    
    % Octaedro
    barplot(a.oct(:, 1), one_or_two(1, -1, i)*a.oct(:, 2)/6, 2.5, clr.oct);
    barplot(a.oct(:, 1), one_or_two(1, -1, i)*a.oct(:, 3)/6, 2.5, clr.oct);
    plot(v.oct1(:, 1), one_or_two(1, -1, i)*v.oct1(:, 2), 'Color', clr.oct, 'LineSmoothing', 'on');
    plot(v.oct2(:, 1), one_or_two(1, -1, i)*v.oct2(:, 2), 'Color', clr.oct, 'LineSmoothing', 'on');
    
    export(pm, sprintf('out/%s_es_bars_oct_%s.dat', name, one_or_two('up', 'dn', i)), [a.oct(:,1), one_or_two(1, -1, i)*a.oct(:,2)/6, one_or_two(1, -1, i)*a.oct(:,3)/6]);
    export(pm, sprintf('out/%s_es_voigt_oct_%s.dat', name, one_or_two('up', 'dn', i)), [v.oct1(:,1), one_or_two(1, -1, i)*v.oct1(:,2), one_or_two(1, -1, i)*v.oct2(:,2)]);
    
    % Total
    plot(v.es_total{i}(:, 1), one_or_two(1, -1, i)*v.es_total{i}(:, 2), '-k', 'LineWidth', 1.5, 'LineSmoothing', 'on');
    export(pm, sprintf('out/%s_es_voigt_%s.dat', name, one_or_two('up', 'dn', i)), [v.es_total{i}(:,1), one_or_two(1, -1, i)*v.es_total{i}(:, 2)]);
end

line([0; 0], [-50, 50], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--');
line([-50; 50], [0, 0], 'LineWidth', 1, 'Color', 'k');

xlim([-15 15]);
ylim([-max(v.es_total{2}(:,2))-0.1 max(v.es_total{1}(:,2))+0.1]);
