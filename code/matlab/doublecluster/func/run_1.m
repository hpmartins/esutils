function o = run_1(is, fs, wc)
    o = cell(3, 20);
    for r = 1:20
        if ~isempty(fs{r})
            [o{1,r}, o{2,r}, o{3,r}] = main_operator(is, fs{r}, r, wc);
        end
    end
end
