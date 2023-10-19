function [y] = gaussian(x0, y0, x, HWHM)
    y = zeros(numel(x), 1);
    for i = 1:numel(x0)
        for f = 1:numel(x)
            y(f) = y0(i)*(sqrt(log(2)/pi)/HWHM)*exp(-log(2)*((x(f)-x0(i))/HWHM)^2);
        end
    end
end
