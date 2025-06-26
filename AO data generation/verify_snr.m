data_hd = csvread('Hd.csv'); 
data_hr = csvread('hr.csv'); 
data_G = csvread('G.csv');
data_theta = csvread('theta.csv'); 
data_W = csvread('W.csv'); 
snr = csvread('snr.csv');

zz = 2;
Hd=reconstructMatrix(1, 2, data_hd(zz,:));
G=reconstructMatrix(25, 2, data_G(zz,:));
Hr=reconstructMatrix(1, 25, data_hr(zz,:)); 

d_theta=reconstructMatrix(25,1, data_theta(zz,:));
d_W = reconstructMatrix(2, 1, data_W(zz,:));

Hd=pd.*Hd_w(:,:,t0);
theta=theta_init(:,:,t0);
Theta=diag(theta');

