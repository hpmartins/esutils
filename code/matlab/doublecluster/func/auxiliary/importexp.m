function exp_data = importexp(file)
    if isempty(file)
        exp_data = [];
    else
        try
            exp_data = importdata(file);
        catch err
            exp_data = [];
        end
    end
end
