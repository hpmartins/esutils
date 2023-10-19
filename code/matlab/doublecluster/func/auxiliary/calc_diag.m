function [nEd1, nEd2, nU1, nU2, nDq1, nDq2, nJ1, nJ2, nEp, npp1, npp2, nEc1, nEc2] = calc_diag(b, info, pm)
    r = info(1, end);

    nEd1 = sum(b(:, 1:10) > 0, 2) - 0*(r > 0)*one_or_two(1, 0, cc(r));
    nEd2 = sum(b(:, 11:20) > 0, 2) - 0*(r > 0)*one_or_two(0, 1, cc(r));
    
    nU1 = nEd1.*(nEd1-1)/2;
    nU2 = nEd2.*(nEd2-1)/2;
    
    nDq1 = 6*sum(b(:, [4 5 9 10]) > 0, 2) - 4*sum(b(:, [1 2 3 6 7 8]) > 0, 2);
    nDq2 = 6*sum(b(:, [14 15 19 20]) > 0, 2) - 4*sum(b(:, [11 12 13 16 17 18]) > 0, 2);
    
    nJ1 = (sum(b(:, 1:5) > 0, 2).*(sum(b(:, 1:5) > 0, 2) - 1))/2 + (sum(b(:, 6:10) > 0, 2).*(sum(b(:, 6:10) > 0, 2) - 1))/2;
    nJ2 = (sum(b(:, 11:15) > 0, 2).*(sum(b(:, 11:15) > 0, 2) - 1))/2 + (sum(b(:, 16:20) > 0, 2).*(sum(b(:, 16:20) > 0, 2) - 1))/2;
    
    nEp = 10 - info(:, 3);

    npp1 = info(:, 4) - info(:, 5);
    npp2 = info(:, 6) - info(:, 7);
    
    nEc1 = sum(b(:, 1:10) == -1, 2) - sum(b(:, 1:10) == 2, 2);
    nEc2 = sum(b(:, 11:20) == -1, 2) - sum(b(:, 11:20) == 2, 2);
end 