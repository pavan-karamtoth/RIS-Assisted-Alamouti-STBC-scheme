clc; clear;
k = 1;
N = 32;
M = 4;
snr_db = 0:5:30;

pe_out_an = zeros(size(snr_db));

% Use high precision to stabilize special functions
digits(50);

for m = 1:numel(snr_db)
    pp = 10^(0.1*snr_db(m));
    lambda = 4*(1+k)/(N^2*pp);
    ebs = lambda / (sin(pi/M)^2);

    sum1 = vpa(0); sum2 = vpa(0);
    Nmax = 80; tol = vpa(1e-14);

    for n = 0:Nmax
        a = -(1+n)/2;
        b = 1 - n/2;
        coef = vpa(k)^n / vpa(factorial(sym(n)));

        % --- term 1 ---
        t1 = coef * (vpa(ebs)^(n/2 + 1)) * exp(-vpa(k) + vpa(ebs)/2) ...
             * whittakerW(a, b, vpa(ebs));
        sum1 = sum1 + t1;

        % --- term 2 ---
        x2 = vpa(3*ebs/4);
        t2 = coef * (x2^(n/2 + 1)) * exp(-vpa(k) + vpa(3*ebs)/8) ...
             * whittakerW(a, b, x2);
        sum2 = sum2 + t2;

        % simple convergence break
        if max(abs([t1, t2])) < tol*max(1, abs(sum1)+abs(sum2))
            break
        end
    end

    pe_out_an(m) = double(sum1/6 + sum2/2);
end

semilogy(snr_db, pe_out_an, 'g--', 'LineWidth', 1.5);
grid on; xlabel('SNR (dB)'); ylabel('SER (approx)');
title('Analytical SER (Whittaker-series), N=32, M=4');
