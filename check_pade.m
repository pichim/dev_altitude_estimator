clc, clear all
%%
Ts = 1/40;
G0 = tf(1,[1/(2*pi*5) 1])
[num, den] = pade(3*Ts, 1);
G = G0 * tf(num, den)

sys0 = ss(G0);
sys  = ss(G);
% sysd = ss(sys.a*Ts + eye(size(sys.a)), sys.b*Ts, sys.c, sys.d, Ts)
sysd = c2d(sys, Ts, 'tustin');
% sysd = c2d(sys, Ts);

figure(1)
step(sys0, sys, sysd), grid on

figure(2)
bode(sys0, sys, sysd), grid on