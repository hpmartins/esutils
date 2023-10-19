so1 = 23.0;
so2 = 6.0;

p1 = 462.6; %posi��o do primeiro pico do espectro 1
p2 = 457.7; %posi��o do primeiro pico do espectro 2


csec1 = 4/7*21.66;
csec2 = 1/2*4.099;
csecp = 6/4*6.977;

N1 = 1;
N2 = 1;
crosssec_cs1 = 1;
crosssec_cs2 = 1;


filename_tag = 0; %numero que identifica o conjunto de parametros usado

wL = 0.5;
aL = 0.5;
wG = 0.7;

x_win = 10;
x_step = 0.05;

p1rs = 1.0;

U1 = 4.0;
D1 = 1.0;

U2 = 4.5;
D2 = 4.0;

kanamori = 0;


Dq1 = 0.37;
Dq2 = 0.25;
J1 = 0.7;
J2 = 0.0;
pps1 = 0.0;
pps2 = 0.0;
pps  = 0.7;
ppp1 = -0.3*pps1;
ppp2 = -0.3*pps2;

Ts1 = 2.2;
Tp1 = - Ts1/2.2;
Ts2 = 4.0;
Tp2 = - Ts2/2.2;


Dc1 = 0.6;
Dc2 = 0.0;
Tc1 = 0.25;
Tc2 = 0.0;

Ep = -6;

Ed1 = D1 + Ep - sum(start(1,1:10),2)*U1;
Ed2 = D2 + Ep - sum(start(1,11:20),2)*U2;
Ec1 = Ed1 - Dc1 + sum(start(1,1:10),2)*U1;
Ec2 = Ed2 - Dc2 + sum(start(1,11:20),2)*U2;

if kanamori
    kanamori_params;
end