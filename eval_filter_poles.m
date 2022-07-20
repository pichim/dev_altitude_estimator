%%

clc, clear variables

n = 3;
w0 = 1;
[b, a] = besself(n, w0);
G = tf(b, a);
poles_G = pole(G);
abs(poles_G)
-cos(atan2(imag(poles_G), real(poles_G)))

s = tf("s");
if n == 2
    abs(roots([1 3 3])) ./ abs(roots(a)) % sqrt(3) for n = 2
    w0b = w0 / sqrt(3);
    H = 3 / (s^2/w0b^2 + 3*s/w0b + 3);
else
    s1 = 0.941600026533207; % 0.149860298638220*2*pi;
    s2 = 1.030544545438434; % 0.164016258482917*2*pi;
    D2 = 0.723540179945206;
    w1 = w0 * s1;
    w2 = w0 * s2;
    H = tf(1, [1/w1 1]) * tf(w2^2, [1 2*D2*w2 w2^2]);
end

opt = bodeoptions;
figure(99)
bode(G, H, logspace(-1, 2, 1e4), opt), grid on

figure(100)
pzmap(G, H), grid on, axis equal

%% 

clc, clear variables

w0 = 2*pi*0.16;
f0 = w0/2/pi;
Ts = 1/1e2;
G = get_filter('pt3', f0, Ts);
poles_G = -1/Ts*log(pole(G));
abs(poles_G) / 2 / pi
cos(atan2(imag(poles_G), real(poles_G)))
s1 = 1.951892355462724; % 0.310653316755175*2*pi;
s2 = 1.953080328553632; % 0.310842388544853*2*pi;
D2 = 0.999999938350424;
w1 = w0 * s1;
w2 = w0 * s2;
H = tf(1, [1/w1 1]) * tf(w2^2, [1 2*D2*w2 w2^2]);

figure(101)
bode(G, H), grid on
