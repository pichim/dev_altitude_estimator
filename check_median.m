clc, clear variables
%%

Ts = 1/40;
N = 5/Ts;
time = (0:N-1).'*Ts;
sig = sin(2*pi*1*time) + 10*randn(N, 1);
sig = cumsum(sig);

n = 3;
sigf_mean = sig;
sigf_median = sig;
for i = 1+n:N
    sigf_mean(i)   = mean(sig(i-n:i));
    sigf_median(i) = median(sig(i-n:i));
end

a = 600*1e-3;
Gf = tf(1-a, [1 -a], Ts)
f_cut = -1/Ts*log(a)/2/pi

sig_f_lp = filter(Gf.num{1}, Gf.den{1}, sigf_median);

figure(1)
stairs(time, [sig, sigf_mean, sigf_median, sig_f_lp]), grid on

[sig(1:10), sigf_mean(1:10), sigf_median(1:10), sig_f_lp(1:10)]

% to check its best to remove noise and not integrate the signal
finddelay(sig, sig_f_lp)