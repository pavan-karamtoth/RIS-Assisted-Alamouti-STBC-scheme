clc;
clear all;

k = 1;
N1 =4; N2=8; N=N1*N2;
M = 4; % M-PSK symbols
snr_th = 100;
var_h1 = 1;
var_h2 = 1;
itr = 1e5;
sigma2 = 1;
Es = 1;
var_n = 1;
mu = 2; var_g = 4;
snr_db = -10:5:20;
pr_out = zeros(1, length(snr_db));

% M-PSK sumbols
M_PSK_angle = ((0:M-1)*2*pi) /M;
%M_PSK_symbols = exp(1j*M_PSK_angle);

for i = 1:length(snr_db)
    % Pr_out = zeros(1, length(snr_db));
    disp(['SNR:', num2str(snr_db(i))])
    p = 10^(0.1*snr_db(i));
    count = 0;

    for j = 1:itr
        g_los = URPA(N1, N2);
        %g_los = mu + var_g*(rand(N,1) + 1j * rand(N,1));
        g = sqrt(k/(1+k)) * g_los + sqrt(1/(1+k)) * wgn(N,1,0,'complex');
      
         a = abs(g);
         theeta = angle(g);
         h1 = sqrt(var_h1) * wgr(N/2, 1, 0, 'complex');
         h2 = sqrt(var_h2) * wgn(N/2, 1, 0, 'complex');
         h = [h1;h2];
        
        % IRS Phase 
        index1 = randi([1 M]);
        index2 = randi([1 M]);
        phai1 = M_PSK_angle(index1);
        phai2 = M_PSK_angle(index2);
        % Modulated Phase angle
        psi1 = theeta(1) + phai1;
        psi2 = theeta(1) + phai2;

        % Modulated Symbols
       z1 = a(1)*sqrt(Es)*exp(1j*psi1);
       z2 = a(1)*sqrt(Es)*exp(1j*psi2);

       % Received equation at Destination in the time slot 1
       y1 = N/2*((h1(1)*z1 + h2(1)*z2)) + sqrt(var_n)*wgn(1,1,0,'complex');
       % Received equation at Destination in the time slot 2
       y2 = N/2*(conj(h2(1))*conj(z1) - conj(h1(1))*conj(z2)) + sqrt(var_n)*wgn(1,1,0,'complex');
       % Alamouti MRC
       h_norm =sqrt(h1(1)^2 + h2(1)^2);
       %z1_hat =(conj(h1(1))/h_norm) *y1 + h2(1)/h_norm * conj(y2);
       %z2_hat = (-h1(1))/h_norm * conj(y2) + conj(h2(1))/h_norm *y1;
       n_tl_1=sqrt(sigma2)*wgn(1,1,0,'complex');
       n_tl_2=sqrt(sigma2)*wgn(1,1,0,'complex');


       z1_hat =N/2 * (h_norm) *z1   + n_tl_1';
       z2_hat =N/2 * (h_norm) *z2   + n_tl_2';
       c1=N/2 * (h_norm) *z1;
       c2=N/2 * (h_norm) *z2;

        % Outage Analysis
        
       snr =  p * ((abs(c1))^2 + (abs(c2))^2)/ sigma2  ;
       %snr= p* ((h_norm)^2* norm(a)^2* N^2) /4;
        

        if snr < snr_th
            count = count + 1;
        end
    end
    pr_out(i) = count/itr;

end

semilogy(snr_db, pr_out, 'r-.+','LineWidth',1);
%figure;
hold on;
%%
% Analytical
pr_out_an = zeros(1,length(snr_db)); 
for m=1:length(snr_db)
    pp = 10^(0.1*snr_db(m));
    Gam = 4*((1+k)/(N^2*pp));
    Sumterm = 0;
    for n = 0:160
    Sumterm =Sumterm + ( (((2* k^n *exp(-k))/(factorial(n))^2 ) *(Gam*snr_th)^((n+1)/2))    *  ((sqrt(Gam*snr_th))*besselk(n,2*sqrt(Gam*snr_th))+besselk(n+1,2*sqrt(Gam*snr_th))));
    end
    pr_out_an(m) = 1-Sumterm;
end
%pr_out_an = 1-Sumterm;
semilogy(snr_db, pr_out_an,'g--','Linewidth', 1);
hold on;





function Reflectinguser = URPA(row, col)
theta_r = 2*pi*rand(1,1);
phai_r = theta_r;
u = cos(theta_r);
v = sin(theta_r) * cos(phai_r);
nx = 0:row-1;
ny = 0:col-1;
Avec1 = exp(1j * pi * nx * u)';
Avec2 = exp(1j * pi * ny * v)';
Reflectinguser = kron(Avec1, Avec2);
end
