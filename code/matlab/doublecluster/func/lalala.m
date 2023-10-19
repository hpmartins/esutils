clf;

x0 = 0;
y0 = 1;
A = [x0; y0]';

wG = 0.1;
wL = 0.1;

x_min = -2.5*(wG+wL);
x_max =  2.5*(wG+wL);
xrange = (x_min:0.01:x_max)';

%% Gaussian
subplot(2, 2, 1); hold on;

title(sprintf('Gaussian HWHM=%.2f', wG));
I  = voigtxr(A, xrange, [1e-10, 0.0, wG]);
yrange = (linspace(min(I), max(I), 100))';

plot(xrange, I, '-k'); 
plot(xrange, max(I)/2);
tmp = xrange(roundn(abs(I-(max(I)/2)), -3) <= 1e-1);
tmp2 = [mean(tmp(tmp<0)), mean(tmp(tmp>0))];
plot(tmp2(1), yrange);
plot(tmp2(2), yrange);
plot(0, yrange);

I_gaussian = gaussian(x0, y0, xrange, wG);
plot(xrange, I_gaussian, '-b');

xlim([min(xrange), max(xrange)]);
ylim([min(yrange), max(yrange)]);
set(gca, 'XTick', [-wG, wG]);

fprintf('Gaussian:\n');
fprintf('HWHM: o que saiu: %f (sigma=%f), o que deveria sair: %f\n', (tmp2(2)-tmp2(1))/2, (tmp2(2)-tmp2(1))/2.35, wG);
fprintf('Altura no máximo: o que saiu: %f, o que deveria sair: %f\n', max(I), y0*(sqrt(log(2)/pi))/wG);
         
%% Lorentzian
subplot(2, 2, 2); hold on;

title(sprintf('Lorentzian HWHM=%.2f', wL));
I  = voigtxr(A, xrange, [wL, 0.0, 1e-10]);
yrange = (linspace(min(I), max(I), 100))';

plot(xrange, I, '-k'); 
plot(xrange, max(I)/2);
tmp = xrange(roundn(abs(I-(max(I)/2)), -3) <= 1e-1);
tmp2 = [mean(tmp(tmp<0)), mean(tmp(tmp>0))];
plot(tmp2(1), yrange);
plot(tmp2(2), yrange);
plot(0, yrange);

I_lorentzian = lorentzian(x0, y0, xrange, wL);
plot(xrange, I_lorentzian+0.1, '-g');

xlim([min(xrange), max(xrange)]);
ylim([min(yrange), max(yrange)]);
set(gca, 'XTick', [-wL, wL]);

fprintf('Lorentzian:\n');
fprintf('HWHM: o que saiu: %f, o que deveria sair: %f\n', (tmp2(2)-tmp2(1))/2, wL);
fprintf('Altura no máximo: o que saiu: %f, o que deveria sair: %f\n', max(I), y0/(3.1415*wL));


%% Voigt
subplot(2, 2, 3); hold on;

title(sprintf('Voigt HWHM_G=%.2f HWHM_L=%.2f', wG, wL));
I  = voigtxr(A, xrange, [wL, 0.0, wG]);
yrange = (linspace(min(I), max(I), 100))';

plot(xrange, I, '-k'); 
plot(xrange, max(I)/2);
tmp = xrange(roundn(abs(I-(max(I)/2)), -3) <= 1e-1);
tmp2 = [mean(tmp(tmp<0)), mean(tmp(tmp>0))];
plot(tmp2(1), yrange);
plot(tmp2(2), yrange);
plot(0, yrange);

I_voigt = voigtpy(x0, y0, xrange, wG, wL);
plot(xrange, I_voigt, '-r');

xlim([min(xrange), max(xrange)]);
ylim([min(yrange), max(yrange)]);
set(gca, 'XTick', [tmp2(1), tmp2(2)]);

fprintf('Voigt:\n');
fprintf('HWHM: o que saiu: %f\n', (tmp2(2)-tmp2(1))/2);
fprintf('Altura no máximo: o que saiu: %f\n', max(I));

%% Tres
subplot(2, 2, 4); hold on;

plot(xrange, I_gaussian, '-b');
plot(xrange, I_lorentzian, '-g');
plot(xrange, I_voigt, '-r');

fprintf('\nPLOT: em preto o nosso voigtxr(), em colorido o que uma gauss/lorentz/voigt deveria ser de fato\n');
fprintf('(se ta sem preto é pq ta em cima\n');