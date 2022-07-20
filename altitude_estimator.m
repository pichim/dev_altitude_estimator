function [x_ret, acc_z_ret] = altitude_estimator(K, nd_pos, wa, acc, pos, quat, Ts)

N = size(acc, 1);
x = zeros(3,1);
x_ret = zeros(N,3);
acc_z_ret = zeros(N,1);

downsample_cntr = 0;
downsample_divider = 1; % if > 1 pos will be downsampled

downsample_pos = 0;

pos_past_idx = 1;
pos_past = zeros(1+nd_pos, 1); % ringbuffer

for i = 1:N
    
    % CEB =
    % [cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)]
    % [cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi)]
    % [        -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta)]
    
    % assume the bias affects only acc z in body frame
    u = acc(i,:).';
    u(3) = u(3) - x(3);
    u = quat2CEB(quat(i,:)) * u;
    u(3) = u(3) - 9.81;
    
    % predict pos, vel and acc bias
    % p = p + v*Ts + 0.5*a*Ts^2;
    % v = v + a*Ts;
    % b = b;
    x = x + ...
        Ts*[x(2) + 0.5*u(3)*Ts; ...
            u(3); ...
            -wa*x(3)];
        
    % update delayed position estimation
    pos_past(pos_past_idx) = x(1);
    pos_past_idx = pos_past_idx + 1;
    if pos_past_idx > 1+nd_pos
        pos_past_idx = 1;
    end
     
    % downample
    if mod(downsample_cntr, downsample_divider) == 0
        downsample_pos = pos(i);
    end
    downsample_cntr = downsample_cntr + 1;
    
    % correction based on delayed position estimate
    x = x + K * Ts * (downsample_pos  - pos_past(pos_past_idx));
    
    x_ret(i,:) = x.';
    acc_z_ret(i,1) = u(3);

end

end

