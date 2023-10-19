clear all
%fofo
% huauhauhha
pm.start = [1,1,1,1,1, 0,0,0,0,0,   0,0,0,0,0, 1,0,0,0,0];
pm.U   = [7.0, 3.0];
pm.Ep  = -6.0;
pm.D   = [2.8, 3.5];
pm.Dq  = [0.15, 0.35];
pm.pds = [3.20, 3.00];
pm.pps = [0.55, 0.55];


pm.pdp = -0.45*pm.pds;
pm.ppp = -0.30*pm.pps;
pm.Ed  = pm.D + pm.Ep;
oct = main_oct(pm);

clf; hold on;

subplot(4, 1, 1);
barplot(oct(:, 1), oct(:, 2), 2, 'b');
set(gca, 'XDir', 'reverse');
ylim([0 max(oct(:, 6))]);

subplot(4, 1, 2);
barplot(oct(:, 1), oct(:, 3), 2, 'r'); 
set(gca, 'XDir', 'reverse');
ylim([0 max(oct(:, 6))]);

subplot(4, 1, 3);
barplot(oct(:, 1), oct(:, 4), 2, 'b');
set(gca, 'XDir', 'reverse');
ylim([0 max(oct(:, 6))]);

subplot(4, 1, 4);
barplot(oct(:, 1), oct(:, 5), 2, 'g'); 
set(gca, 'XDir', 'reverse');
ylim([0 max(oct(:, 6))]);
