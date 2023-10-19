function [nu1, nul1, nu2, nul2] = calc_kanamori(base, nU1, nU2)
    nu1  = sum((base(:, 1:5) + base(:, 6:10)) == 2, 2);
    nul1 = nU1 - nu1;
 
    nu2  = sum((base(:, 11:15) + base(:, 16:20)) == 2, 2);
    nul2 = nU2 - nu2;
end


