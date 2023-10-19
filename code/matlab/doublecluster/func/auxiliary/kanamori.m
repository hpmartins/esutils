function [u1, ul1, j1, u2, ul2, j2] = kanamori(pm)
    F_red  = [pm.F_red1, pm.F_red2];
    F2     = [pm.F21, pm.F22];
    F4     = [pm.F41, pm.F42];

    F2 = F_red.*F2/49;
    F4 = F_red.*F4/441;

    % Calculates the Racah parameters A, B and C.
    C = 35*F4;
    B = F2 - 5*F4;
    A = [pm.U1 pm.U2] + (14/9)*B - (7/9)*C;

    %Calculates the Kanamori parameters u, u' and j
    u = A + 4*B + 3*C;
    ul = A - B + C;
    j = (5/2)*B + C;

    u1 = u(1);
    ul1 = ul(1);
    j1 = j(1);

    u2 = u(2);
    ul2 = ul(2);
    j2 = j(2);
end