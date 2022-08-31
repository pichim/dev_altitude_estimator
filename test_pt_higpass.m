clc, clear all
%%

try
  addpath ../fcn_bib
catch
  addpath fcn_bib_copy
end

Ts = 1/120;

k = 30/100;
f_cut = 1/(2*pi*Ts*(1/k - 1)) / 1.961459176700620;
[G, B, A] = get_filter('pt3', f_cut, Ts);
% f_cut = 1/(2*pi*Ts*(1/k - 1)) / 1.553773974030037
% [G, B, A] = get_filter('pt2', f_cut, Ts);

Gint = tf([Ts 0], [1 -1], Ts);
Ghp = (1 - G);
Gvel2 = (3*Ts*(1 - k))/k * G;
% Gvel2 = (2*Ts*(1 - k))/k * G;

figure(1)
bode(Gint, Ghp, Gint * Ghp, 2*pi*logspace(-2, log10(1/2/Ts), 1e3)), grid on
legend('Gint', '1 - Gpt3', 'Gint * (1 - Gpt3)')

figure(2)
bode(Gint * Ghp, Gvel2, 2*pi*logspace(-2, log10(1/2/Ts), 1e3)), grid on
legend('Gint * (1 - Gpt3)', 'Gpt3 scaled')

figure(3)
step(Gint * Ghp, Gvel2), grid on
legend('Gint * (1 - Gpt3)', 'Gpt3 scaled', 'location', 'best')

%%

% clc, clear all
% 
% syms k z Ts
% 
% % [k 0], [1 (k-1)]
% 
% G = k*z / (z + (k-1))
% Ghp = 1 - G
% Gint = Ts*z / (z - 1)
% Gvel = simplify(Ghp * Gint)
% subs(Gvel, z, 1)
% % -(Ts*(k - 1))/k = Ts*(1 - k)/k

%%

% clc, clear all
% 
% syms k z Ts
% 
% % [k^2 0 0], [1 2*(k-1) (k-1)^2]
% 
% G = k^2*z^2 / (z^2 + 2*(k-1)*z + (k-1)^2)
% Ghp = 1 - G
% Gint = Ts*z / (z - 1)
% Gvel = simplify(Ghp * Gint)
% subs(Gvel, z, 1)
% % -(2*Ts*(k - 1))/k = (2*Ts*(1 - k))/k

%%

% clc, clear all
% 
% syms k z Ts
% 
% % [k^3 0 0 0], [1 3*(k-1) 3*(k-1)^2 (k-1)^3]
% 
% G = k^3*z^3 / (z^3 + 3*(k-1)*z^2 + 3*(k-1)^2*z + (k-1)^3)
% Ghp = 1 - G
% Gint = Ts*z / (z - 1)
% Gvel = simplify(Ghp * Gint)
% subs(Gvel, z, 1)
% % -(3*Ts*(k - 1))/k = (3*Ts*(1 - k))/k