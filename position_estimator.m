function [x_ret, acc_z_ret] = position_estimator(Ts, G, nd_pos, x0, acc, pos, quat)

N = size(acc, 1);
x = x0;
x_ret = zeros(N,3);
acc_z_ret = zeros(N,1);

pos_past_idx = 1;
pos_past = zeros(1+nd_pos, 1); % ringbuffer

% G = [[1, Ts*(1    - k1),   0, Ts^2*(1    - k1), k1]
%      [0,    (1 - Ts*k2),   0, Ts  *(1 - Ts*k2), k2]
%      [0,       - Ts*k3 , a33,       - Ts^2*k3 , k3]];
a33 = G(3,3);

for i = 1:N
    
    % CEB =
    % [cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)]
    % [cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi)]
    % [        -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta)]
    
    % assume the bias affects only acc z in body frame
    u1 = acc(i,:).';
    u1(3) = u1(3) - a33 * x(3);
    u1 = quat2CEB(quat(i,:)) * u1;
    %u1(3) = u1(3) - a33 * x(3);
    u1 = u1(3) - 9.81;
    
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

