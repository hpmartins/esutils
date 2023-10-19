function export(pm, f, x)
    if pm.export
        dlmwrite(f, x, 'delimiter', ' ', 'precision', '%.5f')
    end
end