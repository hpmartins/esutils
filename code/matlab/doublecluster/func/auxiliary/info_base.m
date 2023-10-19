function out = info_base(b)
    b1 = b(:, 1:10); b2 = b(:, 11:end);
    
    ns = min(sum(b == 1, 2)) + any(ismember(-1, b));
    
    d1 = sum(b1 > 0, 2); d2 = sum(b2 > 0, 2);
    C1 = sum(b1 == 2, 2); C2 = sum(b2 == 2, 2);
    c1 = sum(b1 == -1, 2); c2 = sum(b2 == -1, 2);
    L  = d1 + d2 - ns - C1 - C2;
    
    out = [d1 C1 c1 L d2 C2 c2];
end