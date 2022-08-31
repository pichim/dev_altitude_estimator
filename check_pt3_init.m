clc, clear all
%%

try
  addpath ../fcn_bib
catch
  addpath fcn_bib_copy
end

Ts = 1/40;
N = 5/Ts;
time = (0:N-1).'*Ts;
sig = 0.5*randn(N, 1);
sig(100:end) = sig(100:end) - sin(2*pi*1*time(100:end));
sig = cumsum(sig);

% pt3
k = 0.25 * 1.961459176700620;
f_cut = 1/(2*pi*(Ts / k - Ts)) / 1.961459176700620;
[G, B, A] = get_filter('pt3', f_cut, Ts)

f_cut_init = 10*f_cut;
Ninit = 0.8/Ts;

states = sig(1) * eye(3,1);
sig_f = zeros(size(sig));
for i = 1:N
    if i <= Ninit
        f_cut_actual(i) = -(f_cut_init - f_cut)/(Ninit-1)*(i-1) + f_cut_init;
    else
        f_cut_actual(i) = f_cut;
    end
    [~, B, A] = get_filter('pt3', f_cut_actual(i), Ts);
    sig_f(i) = B(1) * sig(i) - A(2:4)*states;
    states = [sig_f(i); states(1:end-1)];
end

figure(1)
stairs(time, f_cut_actual), grid on

figure(2)
stairs(time, [sig, sig_f, filter(B, A, sig)]), grid on
