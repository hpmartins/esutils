function [y] = lorentzian(x0, y0, x, HWHM)
    y = zeros(numel(x), 1);
    for i = 1:numel(x0)
        for f = 1:numel(x)
            y(f) = y0(i)*(HWHM/pi)/((x(f)-x0(i))^2 + HWHM^2);
        end
    end
end
