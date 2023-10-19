init;
%% GS
[b.gs, info.gs, h.gs, e.gs] = run_0(pm, 0, 0);

%% CS
h.cs = cell(1, 2);
for c = 1:2
    h.cs{c} = h.gs;
    h.cs{c} = sparse(h.cs{c} - diag(sum(b.gs(:, anticc(c)) > 0, 2)*pm.U(c)/0.83));
end
e.cs = cell(1, 2);
for c = 1:2
    [e.cs{c}.eva, e.cs{c}.evc] = diagonalize(h.cs{c});
end

%% Transições
a.cs = cell(2, 2);
for c = 1:2
    a.cs{1,c} = [(e.cs{c}.eva - min(e.cs{c}.eva) + pm.cfg.cs_p(c)), pm.weight(c)*pm.csec_core{c}*((e.cs{c}.evc'*e.gs.evc).^2)];
    a.cs{2,c} = a.cs{1,c};
    a.cs{2,c}(:, 1) = a.cs{2,c}(:, 1) + pm.cfg.cs_so(c);
    a.cs{2,c}(:, 2) = a.cs{2,c}(:, 2)/2;
    a.cs{1,c} = a.cs{1,c}(a.cs{1,c}(:, 2)/max(a.cs{1,c}(:, 2)) > 1e-4, :);
    a.cs{2,c} = a.cs{2,c}(a.cs{2,c}(:, 2)/max(a.cs{2,c}(:, 2)) > 1e-4, :);
end

xrange = cell(1, 2);
if (min(a.cs{1,1}(:, 1)) > min(a.cs{1,2}(:, 1)) & min(a.cs{1,1}(:, 1)) < max(a.cs{2,2}(:, 1))) ...
| (min(a.cs{1,2}(:, 1)) > min(a.cs{1,1}(:, 1)) & min(a.cs{1,2}(:, 1)) < max(a.cs{2,1}(:, 1)))
    exp_data = importexp(pm.cfg.cs_exp{1});
    if isempty(exp_data)
        xrange_min = min(min(a.cs{1,1}(:,1)), min(a.cs{1,2}(:,1))) - 10;
        xrange_max = max(max(a.cs{2,1}(:,1)), max(a.cs{2,2}(:,1))) + 10;
    else
        xrange_min = min(exp_data(:, 1)) - 10;
        xrange_max = max(exp_data(:, 1)) + 10;
    end
    xrange{1} = (xrange_min:0.1:xrange_max)';
    xrange{2} = (xrange_min:0.1:xrange_max)';
else
    for c = 1:2
        exp_data = importexp(pm.cfg.cs_exp{c});
        if isempty(exp_data)
            xrange{c} = ((min(a.cs{1,c}(:,1))-10):0.1:(max(a.cs{2,c}(:,1)+10)))';
        else
            xrange{c} = ((min(exp_data(:,1))-10):0.1:(max(exp_data(:,1)+10)))';
        end
    end
end

%% Voigt
v.cs = cell(2, 2);
for c = 1:2
    for p = 1:2
        v.cs{p,c} = voigtxr(a.cs{p,c}(:, 1:2), xrange{c}, [one_or_two(pm.cfg.cs_wL(p), pm.cfg.cs_wL(p+2), c), one_or_two(pm.cfg.cs_aL(p), pm.cfg.cs_aL(p+2), c), one_or_two(pm.cfg.cs_wG(p), pm.cfg.cs_wG(p+2), c)]);
    end
end


%% Plot
if (min(a.cs{1,1}(:, 1)) > min(a.cs{1,2}(:, 1)) & min(a.cs{1,1}(:, 1)) < max(a.cs{2,2}(:, 1))) ...
| (min(a.cs{1,2}(:, 1)) > min(a.cs{1,1}(:, 1)) & min(a.cs{1,2}(:, 1)) < max(a.cs{2,1}(:, 1)))
    clf;
    hold on;

    total_cs = v.cs{1,1} + v.cs{2,1} + v.cs{1,2} + v.cs{2,2};
    total_cs = horzcat(xrange{1}, total_cs);

    % Soma das voigts
    total_max = max(total_cs(:, 2));
    total_cs(:, 2) = normalize(total_cs(:, 2));
    
    % Background & plot
    exp_data = importexp(pm.cfg.cs_exp{1});
    if ~isempty(exp_data) && pm.cfg.cs_bg
        exp_data(:, 2) = normalize(exp_data(:, 2));
        exp_data = sortrows(exp_data, 1);
        scatter(exp_data(:, 1), exp_data(:, 2), '.k');
        if pm.cfg.cs_bg_exp
            exp_bg = [exp_data(:, 1) cumsum(exp_data(:, 2))];
        else
            exp_bg = [total_cs(:, 1) cumsum(total_cs(:, 2))];
        end
        exp_bg(:, 2) =  exp_data(end, 2)*normalize(exp_bg(:, 2));
        exp_bg = interp1(exp_bg(:, 1), exp_bg(:, 2), xrange{1}, 'linear', 'extrap');
        total_cs(:, 2) = normalize(total_cs(:, 2) + exp_bg);
        plot(total_cs(:, 1), total_cs(:, 2), '-k', 'LineSmoothing', 'on');
        [~, idx] = max(total_cs(:, 1));
        plot(total_cs(:, 1), total_cs(idx, 2)*normalize(exp_bg), ':k', 'LineSmoothing', 'on');
        xlim([min(exp_data(:, 1)), max(exp_data(:, 1))]);
    else
        plot(total_cs(:, 1), total_cs(:, 2), '-k', 'LineSmoothing', 'on');
        xlim([min(total_cs(:, 1)), max(total_cs(:, 1))]);
    end
    
    out_voigt = xrange{1};
    for c = 1:2
        a32 = a.cs{1,c};
        a12 = a.cs{2,c};
        out_voigt = padadd(out_voigt, v.cs{1,c});
        out_voigt = padadd(out_voigt, v.cs{2,c});
    end
    if ~isempty(exp_data) && pm.cfg.cs_bg
        out_voigt = padadd(out_voigt, exp_bg);
    end
    out_voigt = padadd(out_voigt, total_cs(:, 2));
    export(pm, sprintf('out/%s_cs_voigt.dat', name), out_voigt);
        
    
    % plotar 3/2, 1/2
    for c = 1:2
        plot(xrange{1}, v.cs{1,c}/(3*total_max), one_or_two('r', 'g', c), 'LineSmoothing', 'on');
        plot(xrange{1}, v.cs{2,c}/(3*total_max), one_or_two('r', 'g', c), 'LineSmoothing', 'on');
        barplot(a.cs{1,c}(:, 1), a.cs{1,c}(:, 2)/3, 2, one_or_two('r', 'g', c));
        barplot(a.cs{2,c}(:, 1), a.cs{2,c}(:, 2)/3, 2, one_or_two('r', 'g', c));
        tmp = padadd(a.cs{1,c}, a.cs{2,c});
        tmp(isnan(tmp)) = 0;
        export(pm, sprintf('out/%s_cs_bars_c%d.dat', name, c), tmp);
    end

    ylim([-0.1, 1.1]);
    set(gca, 'XDir', 'reverse');
    
    title(strcat(name, ' Core Level Spectrum'));
    xlabel('Binding Energy (eV)');
    ylabel('Normalized Intensity');
    set(gca, 'YTickLabel', '');
    set(gca, 'YTick', []);
else
    clf;
    
    for c = 1:2
        subplot(2, 1, c); hold on;
        if isempty(v.cs{c})
            continue;
        end
        
        v_total = v.cs{1,c} + v.cs{2,c};
        total_max = max(v_total);
        v_total = normalize(v_total);
        v.cs{1,c} = v.cs{1,c}/total_max;
        v.cs{2,c} = v.cs{2,c}/total_max;
       
        % Background & plot
        exp_data = importexp(pm.cfg.cs_exp{c});
        if ~isempty(exp_data) && pm.cfg.cs_bg
            exp_data(:, 2) = normalize(exp_data(:, 2));
            exp_data = sortrows(exp_data, 1);
            scatter(exp_data(:, 1), exp_data(:, 2), '.k');
            if pm.cfg.cs_bg_exp
                exp_bg = [exp_data(:, 1) cumsum(exp_data(:, 2))];
            else
                exp_bg = [xrange{c} cumsum(v_total)];
            end
            exp_bg(:, 2) =  exp_data(end, 2)*normalize(exp_bg(:, 2));
            exp_bg = interp1(exp_bg(:, 1), exp_bg(:, 2), xrange{c}, 'linear', 'extrap');
            v_total = normalize(v_total + exp_bg);
            plot(xrange{c}, v_total, '-k', 'LineSmoothing', 'on');
            [~, idx] = max(xrange{c});
            plot(xrange{c}, v_total(idx)*normalize(exp_bg), ':k', 'LineSmoothing', 'on');
            xlim([min(exp_data(:, 1)), max(exp_data(:, 1))]);
            
            export(pm, sprintf('out/%s_cs_total_c%d.dat', name, c), [xrange{c} v_total(idx)*normalize(exp_bg) v_total]);
        else
            plot(xrange{c}, v_total, '-k', 'LineSmoothing', 'on');
            xlim([min(xrange{c}), max(xrange{c})]);
        end
        
        % Plot das barras
        barplot(a.cs{1,c}(:, 1), a.cs{1,c}(:, 2), 2, 'k');
        barplot(a.cs{2,c}(:, 1), a.cs{2,c}(:, 2), 2, 'k');
        export(pm, sprintf('out/%s_cs_bars_c%d.dat', name, c), [a.cs{1,c}; a.cs{2,c}]);
        ylim([-0.1, 1.1]);
        set(gca, 'XDir', 'reverse');
        
        title(strcat(name, ' Core Level Spectrum'));
        xlabel('Binding Energy (eV)');
        ylabel('Normalized Intensity');
        set(gca, 'YTickLabel', '');
        set(gca, 'YTick', []);
    end
end

clear c idx xrange tmp1 tmp2