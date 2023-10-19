function [broadened] = voigtxr(I, x, w)
    if isempty(I)
        broadened = zeros(size(x));
        return;
    end
    
    I = sortrows(I, 1);
    x0 = I(:, 1); y0 = I(:, 2);
    w = num2cell(w); [HWHM_L, a_L, HWHM_G] = w{:}; 
    %HWHM_G = FWHM_G/2;
    %HWHM_L = FWHM_L/2;
    HWHM_L = HWHM_L.*(1+a_L*(x0-x0(1)));
    SIGMA = HWHM_G/sqrt(2*log(2));
    y = zeros(size(x, 1), 1);
    for p = 1:size(x, 1)
        z = ((x(p) - x0) + 1j*HWHM_L)/(sqrt(2)*SIGMA);
        y(p) = y0'*real(Faddeeva_w(z))/(SIGMA*sqrt(2*pi));
    end
    broadened = y;
end

