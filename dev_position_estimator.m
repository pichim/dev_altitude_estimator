clc, clear all
%%

try
  addpath ../fcn_bib
catch
  addpath fcn_bib_copy
end

Ts = 1/120;
f_cut = 0.16;
wa = 2*pi*0.05; % this will be beneficial for horizontal position estimates (x-y)
nd_pos = 9;     % this is somewhat a hack and can lead to a unstable filter

Gint = tf([Ts 0], [1 -1], Ts);

opt = bodeoptions('cstprefs');
opt.YLim = {[1e-7 1e3], [-180 180]};
opt.XLim = [1e-3 1/2/Ts];

% w1c = 2*pi*f_cut
% w2c = 2*pi*f_cut;
% w3c = 2*pi*f_cut;
% p = exp(-[w1c w2c w3c]*Ts);

RC = 1/(2*pi*f_cut);
k  = Ts/(RC + Ts);

a33 = 1/(Ts*wa + 1);

Ad = [[1, Ts,  -Ts^2*a33]
      [0,  1,  -Ts*a33]
      [0,  0,  a33]];
Bd = [Ts^2
      Ts
      0];

% this is the observer gain for 3 real poles (for place you need to make them sligthly different)
w1 = (k-1);
w2 = w1; % * 1.001;
w3 = w1; % * 0.999;
% -> k1 = 1/a33*w1*w2*w3 + 1
% -> k2 = ((2 - k1)*a33 + 1 - w1*w2 - w1*w3 - w2*w3 - k1)/(Ts*a33)
% -> k3 = (- a33^2 + (- w1 - w2 - w3)*a33 - k1 - w1*w2 - w1*w3 - w2*w3 + 1)/(Ts^2*a33^2)
c0 = w1*w2*w3;
c1 = 1 - (w1*w2 + w1*w3 + w2*w3);
k1 = c0/a33 + 1;
k2 = (c1 - k1 + (2 - k1)*a33)/(Ts*a33);
k3 = (c1 - k1 - ((w1 + w2 + w3) + a33)*a33)/(Ts^2*a33^2);
K = [k1; k2; k3]
% K_ = place(Ad.', Ad(1,:).', [-w1, -w2, -w3]).'

% % this is the observer gain for 1 real and 2 compl. conj. poles
% % -> k1 = 1/a33*w1*w2^2 + 1
% % -> k2 = (- w2^2 - 2*D2*w1*w2 - k1 - a33*(k1 - 2) + 1)/(Ts*a33)
% % -> k3 = ((- w1 - 2*D2*w2)*a33 - a33^2 - w2^2 - 2*D2*w1*w2 - k1 + 1)/(Ts^2*a33^2)
% w1c = 2*pi*f_cut;
% w2c = 2*pi*f_cut;
% D2c = 1/sqrt(2);  % D = 1/2/Q
% p = exp([-w1c, -w2c*D2c + 1i*w2c*sqrt(1 - D2c^2), -w2c*D2c - 1i*w2c*sqrt(1 - D2c^2)]*Ts);
% w1 = -abs(p(1));
% w2 = -abs(p(2));
% D2 = cos(angle(p(2)));
% k1 = 1/a33*w1*w2^2 + 1;
% k2 = (1 - k1 - w2*(w2 + 2*D2*w1) + (2 - k1)*a33)/(Ts*a33);
% k3 = (1 - k1 - w2*(w2 + 2*D2*w1) - ((w1 + 2*D2*w2)*a33 + a33)*a33)/(Ts^2*a33^2);
% K = [k1; k2; k3]
% K_ = place(Ad.', Ad(1,:).', [-w1, -w2*D2 + 1i*w2*sqrt(1 - D2^2), -w2*D2 - 1i*w2*sqrt(1 - D2^2)]).'


%  x1_k  + Ts*(1    - k1) * x2_k              + Ts^2*(1    - k1) * (u1_k - a33*x3_k) + k1 * (u2_k - x1_k)
%             (1 - Ts*k2) * x2_k              + Ts  *(1 - Ts*k2) * (u1_k - a33*x3_k) + k2 * (u2_k - x1_k)
%                - Ts*k3  * x2_k + a33 * x3_k -         Ts^2*k3  * (u1_k - a33*x3_k) + k3 * (u2_k - x1_k)
% G = [[1, Ts*(1    - k1),   0, Ts^2*(1    - k1), k1]
%      [0,    (1 - Ts*k2),   0, Ts  *(1 - Ts*k2), k2]
%      [0,       - Ts*k3 , a33,       - Ts^2*k3 , k3]];
a12 = Ts*(1 - k1);
a22 = (1 - Ts*k2);
a32 = -Ts*k3;
G = [[1, a12,   0, Ts*a12, k1]
     [0, a22,   0, Ts*a22, k2]
     [0, a32, a33, Ts*a32, k3]];

sys_lin = ss(Ad - K*Ad(1,:), [Bd - K*Bd(1,:), K], ...
             Ad - K*Ad(1,:), [Bd - K*Bd(1,:), K], Ts);

% this is the same like above only if nd_pos = 0
[aa, bb, cc, dd] = dlinmod('position_estimator_G');
sys = ss(aa, bb, cc, dd, Ts);

% they are the complementary parts (sum up to 1)
G_inp = sys(1,1);
G_out = sys(1,2);
G_hp  = G_inp / Gint / Gint;
G_lp  = G_out;

G_dinp = sys(2,1);
G_dout = sys(2,2);
G_dhp = G_dinp;
G_dlp = G_dout * Gint * Gint;

G_ddinp = sys(3,1);
G_ddout = sys(3,2);
G_ddhp = G_ddinp;
G_ddlp = G_ddout * Gint * Gint;

figure(11)
bode(G_hp, G_lp, G_hp + G_lp, opt), grid on

figure(22)
bode(G_dhp, G_dlp, G_dhp + G_dlp, opt), grid on

figure(33)
bode(G_ddhp, G_ddlp, G_ddhp + G_ddlp, opt), grid on

figure(44)
bode(sys(1,1), sys(1,2), sys(2,1), sys(2,2), sys(3,1), sys(3,2), opt), grid on

figure(55)
subplot(121)
step(sys(1,1), sys(1,2), 10), grid on
subplot(122)
step(sys(2,1), 'r', sys(2,2), 'c', 10), grid on

damp(sys(1,1))

%%

% clc, clear all
% syms Ts a33 k1 k2 k3 x1_k x2_k x3_k u1_k u2_k
% 
% Ad = [[1, Ts,  -Ts^2*a33]
%       [0,  1,  -Ts*a33]
%       [0,  0,  a33]];
% Bd = [Ts^2
%       Ts
%       0];
% K = [k1; k2; k3]
% 
% x_k = [x1_k; x2_k; x3_k]
% u_k = [u1_k; u2_k]
% 
% x_kplus1 = (Ad - K*Ad(1,:)) * x_k + [Bd - K*Bd(1,:), K] * u_k
% 
% % + k1*u2_k - x3_k*(Ts^2*a33 - Ts^2*a33*k1) - x1_k*(k1 - 1) + x2_k*(Ts - Ts*k1) - u1_k*(Ts^2*k1 - Ts^2)
% %        + u1_k*(- k2*Ts^2 + Ts) - x3_k*(- a33*k2*Ts^2 + a33*Ts) + k2*u2_k - k2*x1_k - x2_k*(Ts*k2 - 1)
% %                            + k3*u2_k - k3*x1_k + x3_k*(a33*k3*Ts^2 + a33) - Ts*k3*x2_k - Ts^2*k3*u1_k
% 
% %  (1 - k1) * x1_k + Ts*(1    - k1) * x2_k            + Ts^2*(1    - k1) * (u1_k - a33*x3_k) + k1*u2_k
% %      -k2  * x1_k +    (1 - Ts*k2) * x2_k            + Ts  *(1 - Ts*k2) * (u1_k - a33*x3_k) + k2*u2_k
% %      -k3  * x1_k         - Ts*k3  * x2_k + a33*x3_k -         Ts^2*k3  * (u1_k - a33*x3_k) + k3*u2_k
% 
% %  x1_k  + Ts*(1    - k1) * x2_k              + Ts^2*(1    - k1) * (u1_k - a33*x3_k) + k1 * (u2_k - x1_k)
% %             (1 - Ts*k2) * x2_k              + Ts  *(1 - Ts*k2) * (u1_k - a33*x3_k) + k2 * (u2_k - x1_k)
% %                - Ts*k3  * x2_k + a33 * x3_k -         Ts^2*k3  * (u1_k - a33*x3_k) + k3 * (u2_k - x1_k)
% 
% simplify(x1_k  + Ts*(1    - k1) * x2_k              + Ts^2*(1    - k1) * (u1_k - a33*x3_k) + k1 * (u2_k - x1_k) - x_kplus1(1))
% simplify(           (1 - Ts*k2) * x2_k              + Ts  *(1 - Ts*k2) * (u1_k - a33*x3_k) + k2 * (u2_k - x1_k) - x_kplus1(2))
% simplify(              - Ts*k3  * x2_k + a33 * x3_k -         Ts^2*k3  * (u1_k - a33*x3_k) + k3 * (u2_k - x1_k) - x_kplus1(3))
% 
% % collect(x_kplus1, [x1_k x2_k x3_k u1_k u2_k])
% 
% % k1*u2_k - x3_k*(Ts^2*a33 - Ts^2*a33*k1) - x1_k*(k1 - 1) + x2_k*(Ts - Ts*k1) - u1_k*(Ts^2*k1 - Ts^2)
% %        u1_k*(- k2*Ts^2 + Ts) - x3_k*(- a33*k2*Ts^2 + a33*Ts) + k2*u2_k - k2*x1_k - x2_k*(Ts*k2 - 1)
% %                            k3*u2_k - k3*x1_k + x3_k*(a33*k3*Ts^2 + a33) - Ts*k3*x2_k - Ts^2*k3*u1_k

%%

% clc, clear all
% syms wa Ts s k1 k2 k3 z w1 w2 w3
% 
% A = [[0 1 0]; [0 0 -1]; [0 0 -wa]];
% B = [0 1 0].';
% Ad = (eye(3) - Ts*A)^-1
% Bd = Ad * B*Ts
% Cd = Ad
% Dd = Bd
% 
% % Ad =
% % [1, Ts, -Ts^2/(Ts*wa + 1)]
% % [0,  1,   -Ts/(Ts*wa + 1)]
% % [0,  0,     1/(Ts*wa + 1)]
% %  
% % Bd =
% % Ts^2
% %   Ts
% %    0
% %  
% % Cd =
% % [1, Ts, -Ts^2/(Ts*wa + 1)]
% % [0,  1,   -Ts/(Ts*wa + 1)]
% % [0,  0,     1/(Ts*wa + 1)]
% %  
% % Dd =
% % Ts^2
% %   Ts
% %    0
% 
% K = [k1; k2; k3]
% 
% collect(det(z*eye(3) - (Ad - K*Cd(1,:))), 'z')
% % -> z^3 + ((k1 + Ts*k2 - 2*Ts*wa - Ts^2*k3 + Ts*k1*wa + Ts^2*k2*wa - 3)*z^2)/(Ts*wa + 1) - ((2*k1 + Ts*k2 - Ts*wa + Ts*k1*wa - 3)*z)/(Ts*wa + 1) + (k1 - 1)/(Ts*wa + 1)
% 
% collect((w1*w2*w3)*(1/w1*z + 1)*(1/w2*z + 1)*(1/w3*z + 1), 'z')
% % -> z^3 + (w1 + w2 + w3)*z^2 + (w1*w2 + w1*w3 + w2*w3)*z + w1*w2*w3
% 
% % (k1 + Ts*k2 - 2*Ts*wa - Ts^2*k3 + Ts*k1*wa + Ts^2*k2*wa - 3)/(Ts*wa + 1)  != w1 + w2 + w3
% % -(2*k1 + Ts*k2 - Ts*wa + Ts*k1*wa - 3)/(Ts*wa + 1)                        != w1*w2 + w1*w3 + w2*w3
% % (k1 - 1)/(Ts*wa + 1)                                                      != w1*w2*w3
% % -> k1 = w1*w2*w3*(Ts*wa + 1) + 1
% 
% eqns = [(k1 + Ts*k2 - 2*Ts*wa - Ts^2*k3 + Ts*k1*wa + Ts^2*k2*wa - 3)/(Ts*wa + 1)  == w1 + w2 + w3, ...
%         -(2*k1 + Ts*k2 - Ts*wa + Ts*k1*wa - 3)/(Ts*wa + 1)                        == w1*w2 + w1*w3 + w2*w3];
% S = solve(eqns,[k2 k3])
% collect(S.k2, [wa Ts w1])
% % -> k2 = ((- w2 - w3)*wa*Ts*w1 + (1 - w2*w3 - k1)*wa*Ts + (- w2 - w3)*w1 + 3 - w2*w3 - 2*k1)/Ts
% collect(S.k3, [wa Ts w1])
% % -> k3 = ((- w2 - w3)*wa^2*Ts^2*w1 + (1 - w2*w3 - k1)*wa^2*Ts^2 + (- 2*w2 - 2*w3 - 1)*wa*Ts*w1 + (2 - w2 - w3 - 2*w2*w3 - 2*k1)*wa*Ts + (- w2 - w3 - 1)*w1 - k1 - w2 - w3 - w2*w3)/Ts^2
% 
% 
% 
% syms a12 a13 a23 a33
% 
% Ad = [[1, a12, a13]
%       [0,   1, a23]
%       [0,   0, a33]];
% Cd = Ad;
% 
% collect(det(z*eye(3) - (Ad - K*Cd(1,:))), 'z')
% % -> z^3 + (k1 - a33 + a12*k2 + a13*k3 - 2)*z^2 + (2*a33 - k1 - a13*k3 - a33*k1 + a12*a23*k3 - a12*a33*k2 + 1)*z - a33 + a33*k1
% 
% collect((w1*w2*w3)*(1/w1*z + 1)*(1/w2*z + 1)*(1/w3*z + 1), 'z')
% % -> z^3 + (w1 + w2 + w3)*z^2 + (w1*w2 + w1*w3 + w2*w3)*z + w1*w2*w3
% 
% % (k1 - a33 + a12*k2 + a13*k3 - 2)                             != w1 + w2 + w3
% % (2*a33 - k1 - a13*k3 - a33*k1 + a12*a23*k3 - a12*a33*k2 + 1) != w1*w2 + w1*w3 + w2*w3
% % - a33 + a33*k1                                               != w1*w2*w3
% % -> k1 = 1/a33*w1*w2*w3 + 1
% 
% eqns = [(k1 - a33 + a12*k2 + a13*k3 - 2)                             == w1 + w2 + w3, ...
%         (2*a33 - k1 - a13*k3 - a33*k1 + a12*a23*k3 - a12*a33*k2 + 1) == w1*w2 + w1*w3 + w2*w3];
% S = solve(eqns,[k2 k3])
% collect(S.k2, [a12 a13 a23 a33 w1 w2 w3 k1])
% % -> k2 = (a12*a23*a33 + a12*a23*w1 + a12*a23*w2 + a12*a23*w3 - a12*a23*k1 + 2*a12*a23 - a13*a33*k1 + a13*a33 - a13*w1*w2 - a13*w1*w3 - a13*w1 - a13*w2*w3 - a13*w2 - a13*w3 - a13)/(a12^2*a23 + a12*a13*a33 - a12*a13)
% collect(S.k3, [a12 a13 a23 a33 w1 w2 w3 k1])
% % -> k3 = (a33^2 + a33*w1 + a33*w2 + a33*w3 + w1*w2 + w1*w3 + w2*w3 + k1 - 1)/(a12*a23 + a13*a33 - a13)
% 
% 
% 
% % a33 = 1/(Ts*wa + 1);
% % a12 = Ts;
% % a13 = -Ts^2*a33;
% % a23 = -Ts*a33;
% Ad = [[1, Ts, -Ts^2*a33]
%       [0,  1,   -Ts*a33]
%       [0,  0,       a33]];
% Cd = Ad;
% 
% collect(det(z*eye(3) - (Ad - K*Cd(1,:))), 'z')
% % -> z^3 + (- a33*k3*Ts^2 + k2*Ts - a33 + k1 - 2)*z^2 + (2*a33 - k1 - a33*k1 - Ts*a33*k2 + 1)*z - a33 + a33*k1
% 
% collect((w1*w2*w3)*(1/w1*z + 1)*(1/w2*z + 1)*(1/w3*z + 1), 'z')
% % -> z^3 + (w1 + w2 + w3)*z^2 + (w1*w2 + w1*w3 + w2*w3)*z + w1*w2*w3
% 
% % (- a33*k3*Ts^2 + k2*Ts - a33 + k1 - 2) != w1 + w2 + w3
% % (2*a33 - k1 - a33*k1 - Ts*a33*k2 + 1)  != w1*w2 + w1*w3 + w2*w3
% % - a33 + a33*k1                         != w1*w2*w3
% % -> k1 = 1/a33*w1*w2*w3 + 1
% 
% eqns = [(- a33*k3*Ts^2 + k2*Ts - a33 + k1 - 2) == w1 + w2 + w3, ...
%         (2*a33 - k1 - a33*k1 - Ts*a33*k2 + 1) == w1*w2 + w1*w3 + w2*w3];
% S = solve(eqns,[k2 k3])
% collect(S.k2, a33)
% % -> k2 = ((2 - k1)*a33 + 1 - w1*w2 - w1*w3 - w2*w3 - k1)/(Ts*a33)
% collect(S.k3, a33)
% % -> k3 = (- a33^2 + (- w1 - w2 - w3)*a33 - k1 - w1*w2 - w1*w3 - w2*w3 + 1)/(Ts^2*a33^2)
% 
% 
% syms D2
% % a33 = 1/(Ts*wa + 1);
% % a12 = Ts;
% % a13 = -Ts^2*a33;
% % a23 = -Ts*a33;
% Ad = [[1, Ts, -Ts^2*a33]
%       [0,  1,   -Ts*a33]
%       [0,  0,       a33]];
% Cd = Ad;
% 
% collect(det(z*eye(3) - (Ad - K*Cd(1,:))), 'z')
% % -> z^3 + (- a33*k3*Ts^2 + k2*Ts - a33 + k1 - 2)*z^2 + (2*a33 - k1 - a33*k1 - Ts*a33*k2 + 1)*z - a33 + a33*k1
% 
% collect(w1*(1/w1*z + 1)*(z^2 + 2*D2*w2*z + w2^2), 'z')
% % -> z^3 + (w1 + 2*D2*w2)*z^2 + (w2^2 + 2*D2*w1*w2)*z + w1*w2^2
% 
% % (- a33*k3*Ts^2 + k2*Ts - a33 + k1 - 2) != w1 + 2*D2*w2
% % (2*a33 - k1 - a33*k1 - Ts*a33*k2 + 1)  != w2^2 + 2*D2*w1*w2
% % - a33 + a33*k1                         != w1*w2^2
% % -> k1 = 1/a33*w1*w2^2 + 1
% 
% eqns = [(- a33*k3*Ts^2 + k2*Ts - a33 + k1 - 2) == w1 + 2*D2*w2, ...
%         (2*a33 - k1 - a33*k1 - Ts*a33*k2 + 1) == w2^2 + 2*D2*w1*w2];
% S = solve(eqns,[k2 k3])
% collect(S.k2, a33)
% % -> k2 = (- w2^2 - 2*D2*w1*w2 - k1 - a33*(k1 - 2) + 1)/(Ts*a33)
% collect(S.k3, a33)
% % -> k3 = ((- w1 - 2*D2*w2)*a33 - a33^2 - w2^2 - 2*D2*w1*w2 - k1 + 1)/(Ts^2*a33^2)



% T = Cd; % xn = T * x <-> x = T^-1 * xn
% Tinv = simplify(Cd^-1)
% Ad = simplify(T * Ad * T^-1)
% Bd = simplify(T * Bd)
% Cd = simplify(Cd * T^-1)
% Dd = simplify(Dd)