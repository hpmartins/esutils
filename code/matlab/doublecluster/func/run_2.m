function a = run_2(is, o, fs)
    a = cell(3, 20);
    for r = 1:20
        for p = 1:3
            if ~isempty(o{p,r})
                a{p,r} = [(fs{r}.eva - is.eva) ((is.evc' * full(o{p,r}) * fs{r}.evc).^2)' cc(r)*ones(size(fs{r}.eva)) r*ones(size(fs{r}.eva))];
            end
        end
    end
end