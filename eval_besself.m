clc, clear variables

w0 = 1;
[b, a] = besself(2, w0);
G = tf(b, a);

s = tf("s");
abs(roots([1 3 3])) ./ abs(roots(a)) % sqrt(3) for n = 2
w0b = w0 / sqrt(3);
H = 3 / (s^2/w0b^2 + 3*s/w0b + 3);

opt = bodeoptions;
figure(99)
bode(G, H, logspace(-1, 2, 1e4), opt), grid on

figure(100)
pzmap(G, H), grid on, axis equal