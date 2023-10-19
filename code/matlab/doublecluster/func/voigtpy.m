function [y] = voigtpy(x0, y0, x, HWHM_G, HWHM_L)
    SIGMA = HWHM_G/sqrt(2*log(2));
    y = zeros(numel(x), 1);
    for i = 1:numel(x0)
        for f = 1:numel(x)
            z = ((x(f)-x0(i)) + 1j*HWHM_L)/SIGMA/sqrt(2);
            y(f) = y0(i)*real(Faddeeva_w(z))/SIGMA/sqrt(2*pi);
        end
    end
end
