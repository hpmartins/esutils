clear all

%%% Parameters
file = 'NdCrTiO5-AFMup.dat';
w = 0.30; % Gaussian
%%%

data = importdata(file);
x0 = data(:, 1); y0 = data(:, 6);
x = (-10:0.1:30)';
y = zeros(size(x, 1), 1);

for q = 1:size(y0, 1)
    y = y + y0(q)*exp(-(x-x0(q)).^2/(2*w^2));
end
y = (y/trapz(x,y))*trapz(x0,y0);

clf
hold on;
plot(x0, y0, '-r');
plot(x, y, '-g');

[~, f, ~] = fileparts(file);
dlmwrite(strcat(f, '_broadened.dat'), [x y], ' ');