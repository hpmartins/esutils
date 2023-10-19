function out = slater(pm, r, c, t1, t2, o1, o2)
    [l, m, n] = split(r/norm(r));
    
    if (t1 == 1 && t2 == 1)
        p_p;
    elseif (t1 == 1 && t2 == 2)
        p_d(o1, o2);
    elseif (t1 == 2 && t2 == 1)
        p_d(o2, o1);
    else
        out = 0;
    end

    function p_p
        lmn = [l m n];
        if o1 == o2
            out = (lmn(o1)^2)*pm.pps(c) + (1 - (lmn(o1))^2)*pm.ppp(c);
        else
            out = lmn(o1)*lmn(o2)*(pm.pps(c)-pm.ppp(c));
        end
    end

    function p_d(p, d)
        if (p == 1)
          if (d == 1)
             out = l*(n^2-(1/2)*(l^2+m^2))*pm.pds(c)-sqrt(3)*l*(n^2)*pm.pdp(c);
          elseif (d == 2)
             out = (sqrt(3)/2)*l*(l^2-m^2)*pm.pds(c)+l*(1-l^2+m^2)*pm.pdp(c);
          elseif (d == 3)
             out = sqrt(3)*(l^2)*m*pm.pds(c)+m*(1-2*l^2)*pm.pdp(c);
          elseif (d == 4)
             out = sqrt(3)*l*m*n*pm.pds(c)-2*l*m*n*pm.pdp(c);
          elseif (d == 5)
             out = sqrt(3)*(l^2)*n*pm.pds(c)+n*(1-2*l^2)*pm.pdp(c);
          end
        elseif (p == 2)
          if (d == 1)
            out = m*(n^2-(1/2)*(l^2+m^2))*pm.pds(c)-sqrt(3)*m*(n^2)*pm.pdp(c);
          elseif (d == 2)
            out = (sqrt(3)/2)*m*(l^2-m^2)*pm.pds(c)-m*(1+l^2-m^2)*pm.pdp(c);
          elseif (d == 3)
            out = sqrt(3)*(m^2)*l*pm.pds(c)+l*(1-2*m^2)*pm.pdp(c);
          elseif (d == 4)
            out = sqrt(3)*(m^2)*n*pm.pds(c)+n*(1-2*m^2)*pm.pdp(c);
          elseif (d == 5)
            out = sqrt(3)*l*m*n*pm.pds(c)-2*l*m*n*pm.pdp(c);
          end
          elseif (p == 3)
          if (d == 1)
            out = n*(n^2-(1/2)*(l^2+m^2))*pm.pds(c)+sqrt(3)*n*(l^2+m^2)*pm.pdp(c);
          elseif (d == 2)
            out = (sqrt(3)/2)*n*(l^2-m^2)*pm.pds(c)-n*(l^2-m^2)*pm.pdp(c);
          elseif (d == 3)
            out = sqrt(3)*l*m*n*pm.pds(c)-2*l*m*n*pm.pdp(c);
          elseif (d == 4)
            out = sqrt(3)*(n^2)*m*pm.pds(c)+m*(1-2*n^2)*pm.pdp(c);
          elseif (d == 5)
            out = sqrt(3)*(n^2)*l*pm.pds(c)+l*(1-2*n^2)*pm.pdp(c);
          end
        end
    end
end

function varargout = split(a)
    for k = 1:nargout
      varargout{k} = a(k);
    end
end