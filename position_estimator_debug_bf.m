function [x_ret, acc_z_ret] = position_estimator_debug_bf(f_cut, f_a, nd_pos, Ts, acc, pos)

N = size(acc, 1);
x = zeros(3,1);
x_ret = zeros(N,3);
acc_z_ret = zeros(N,1);

pos_past_idx = 1;
pos_past = zeros(1+nd_pos, 1); % ringbuffer

% define location of discrete poles
k  = Ts / (1 / (2 * pi * f_cut) + Ts);
w1 = k - 1;
w2 = w1;
w3 = w1;

wa = 2 * pi * f_a;
a33 = 1 / (Ts*wa + 1);

% calculate observer gain
c0 = w1 * w2 * w3;
c1 = 1 - (w1 * w2 + w1 * w3 + w2 * w3);
c2 = Ts * a33;
k1 = c0 / a33 + 1;
k2 = (c1 - k1 + (2 - k1) * a33) / c2;
k3 = (c1 - k1 - ((w1 + w2 + w3) + a33) * a33) / (c2*c2);

a12 = Ts * (1 - k1);
a22 = (1 - Ts * k2);
a32 = -Ts * k3;

G = [[1, a12,   0, Ts*a12, k1]
     [0, a22,   0, Ts*a22, k2]
     [0, a32, a33, Ts*a32, k3]];
 
for i = 1:N
    
    u1 = acc(i) - a33 * x(3) - 981.0;
    
    % update delayed position estimation
    pos_past(pos_past_idx) = x(1);
    pos_past_idx = pos_past_idx + 1;
    if pos_past_idx > 1+nd_pos
        pos_past_idx = 1;
    end
    u2 = pos(i) - pos_past(pos_past_idx);
    
    %  x1_k  + Ts*(1    - k1) * x2_k              + Ts^2*(1    - k1) * (u1_k - a33*x3_k) + k1 * (u2_k - x1_k)
    %             (1 - Ts*k2) * x2_k              + Ts  *(1 - Ts*k2) * (u1_k - a33*x3_k) + k2 * (u2_k - x1_k)
    %                - Ts*k3  * x2_k + a33 * x3_k -         Ts^2*k3  * (u1_k - a33*x3_k) + k3 * (u2_k - x1_k)
    y = G * [x; u1; u2];
    x = y(1:3);
    
    x_ret(i,:) = x.';
    acc_z_ret(i,1) = u1;

end

end

