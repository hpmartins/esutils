% Init
clear all

dir = 'SFO113n';
ch  = '1';

%%%%%%%%%%%%%%%%%%%%%%
% Parametros
V = 0.9;
d = 0.5; % k-k

V_t2g = -0.13*V;
V_eg  = -0.36*V;
V_d   = -0.00*V; 
V_sp  = -0.00*V;

desloc = 528.1;
%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
d_t2g = importdata(strcat(dir, '/DATA-t2g.dat'));
d_eg = importdata(strcat(dir, '/DATA-eg.dat'));
d_sp = importdata(strcat(dir, '/DATA-sp.dat'));
d_d = importdata(strcat(dir, '/DATA-d.dat'));
XAS = importdata(strcat(dir, '/DATA-XAS.dat'));
%%%%%%%%%%%%%%%%%%%%%%
E = d_t2g(:,1) + desloc;

%%%%%%%%%%%%%%%%%%%%%%
% C�lculo de G0
% G0: G0...
%%%%%%%%%%%%%%%%%%%%%%
G0_t2g = zeros(length(E), 1);
G0_eg  = zeros(length(E), 1);
G0_d   = zeros(length(E), 1);
G0_sp  = zeros(length(E), 1);

for p = 1:length(E)
  x = d_t2g(p,1);
  G0_t2g(p,1) = sum(d_t2g(:,2)./((x-d_t2g(:,1))+1i*d));
  G0_eg(p,1) = sum(d_eg(:,2)./((x-d_eg(:,1))+1i*d));
  G0_d(p,1) = sum(d_d(:,2)./((x-d_d(:,1))+1i*d));
  G0_sp(p,1) = sum(d_sp(:,2)./((x-d_sp(:,1))+1i*d));
end
%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
% C�lculo da nova DOS
% A: nova DOS
%%%%%%%%%%%%%%%%%%%%%%
A_t2g = (-1/pi)*imag((G0_t2g(:,1)./(1-G0_t2g(:,1)*V_t2g)));
A_eg  = (-1/pi)*imag((G0_eg(:,1)./(1-G0_eg(:,1)*V_eg)));
A_d   = (-1/pi)*imag((G0_d(:,1)./(1-G0_d(:,1)*V_d)));
A_sp  = (-1/pi)*imag((G0_sp(:,1)./(1-G0_sp(:,1)*V_sp)));
Atmp  = A_t2g + A_eg + A_sp + A_d;
%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
% Gaussiana
gauss_w = 0.5; % gauss
alfa    = 0;
%%%%%%%%%%%%%%%%%%%%%%
gauss_n = length(Atmp);
gauss_x_min = round(E(1));
gauss_x_max = round(E(gauss_n));
gauss_out_x = gauss_x_min:0.1:gauss_x_max;
gauss_out_y = zeros(1,length(gauss_out_x));

for q = 1:gauss_n
    gauss_w2 = gauss_w*(1+alfa*(Atmp(q)-Atmp(1)));
    gauss_out_y = gauss_out_y + (Atmp(q)/(gauss_w2*sqrt(2*pi)))*gaussmf(gauss_out_x, [gauss_w2 E(q)]);
end
gauss_out_y = (gauss_out_y-min(gauss_out_y))/max(gauss_out_y);
%%%%%%%%%%%%%%%%%%%%%%



%
%plot(XAS(:,1),XAS(:,2),'--',dt(:,1),A)
%
%   E,A_t2g/3.5,'-b', E,A_eg_up/3.5,'-g', E,A_eg_dn/3.5,'-b', E,A_d/2.2,'-r', E,A_sp/2.5,'-g'
plot(XAS(:,1),XAS(:,2),'.k',  gauss_out_x,gauss_out_y,'-k', E,A_t2g/2,'-b', E,A_eg/2.5,'-g', E,A_d/3,'-c', E,A_sp/3,'-m')
axis([520 555 0 1])
out(1:length(gauss_out_x),1) = transpose(gauss_out_x);
out(1:length(gauss_out_y),2) = transpose(gauss_out_y);
dlmwrite(strcat('DATA-', dir, '-', ch, '.dat'), out, 'delimiter', ' ');