%% Includes
addpath('./yaml');
addpath('./func');
addpath('./func/auxiliary');
addpath('./func/utils');
addpath('./func/faddeeva');

%% Input
%clearvars -except name
pm = ReadYaml(strcat(name, '.in'), 0, 1);
pm = set_default_values(pm);

%% Model parameters
pm.pps = [pm.pps1, pm.pps2];
pm.ppp = -0.3*pm.pps;

pm.Ts  = [pm.Ts1, pm.Ts2];
pm.pds = pm.Ts/sqrt(3);
pm.pdp = -0.45*pm.pds;
pm.Tp  = 2*pm.pdp; 

pm.Ed1 = pm.D1 + pm.Ep - sum(pm.start(1,1:10),2)*pm.U1;
pm.Ed2 = pm.D2 + pm.Ep - sum(pm.start(1,11:20),2)*pm.U2;
pm.Ec1 = pm.Ed1 - pm.Dc1 + sum(pm.start(1,1:10),2)*pm.U1;
pm.Ec2 = pm.Ed2 - pm.Dc2 + sum(pm.start(1,11:20),2)*pm.U2;

pm.D = [pm.D1, pm.D2];
pm.U = [pm.U1, pm.U2];
pm.Ed = [pm.Ed1, pm.Ed2];
pm.Ec = [pm.Ec1, pm.Ec2];
pm.Dq = [pm.Dq1, pm.Dq2];

%% Kanamori
[pm.u1, pm.ul1, pm.j1, pm.u2, pm.ul2, pm.j2] = kanamori(pm);
if (pm.J1 < 0)
    pm.J1 = pm.j1;
end
if (pm.J2 < 0)
    pm.J2 = pm.j2;
end

%% Cross-sections and weight
pm.weight = pm.weight/max(pm.weight);
tmpcsec = max([pm.cfg.vb_csec1, pm.cfg.vb_csec2, pm.cfg.vb_csecp]);
pm.cfg.vb_csec1 = pm.cfg.vb_csec1/tmpcsec;
pm.cfg.vb_csec2 = pm.cfg.vb_csec2/tmpcsec;
pm.cfg.vb_csecp = pm.cfg.vb_csecp/tmpcsec;
pm.csec = {pm.cfg.vb_csec1, pm.cfg.vb_csec2; pm.cfg.vb_csecp, pm.cfg.vb_csecp};
pm.csec_core = {pm.cfg.core_csec1 , pm.cfg.core_csec2}; %adicionei csec no core