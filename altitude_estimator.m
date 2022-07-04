function [x_ret, acc_z_ret] = altitude_estimator(K, wa, acc, pos, quat, Ts)

N = size(acc, 1);
x = zeros(3,1);
x_ret = zeros(N,3);
acc_z_ret = zeros(N,2);

downsample_cntr = 0;
downsample_divider = 10;

downsample_pos = 0;

% p = p + v*Ts + 0.5*a*Ts^2;
% v = v + a*Ts;
% b = b;
for i = 1:N
    
%     % directly take acc
%     u = acc(i,:).';
%     u(3) = u(3) - 9.81 - x(3);

    % assume the bias affects only acc z in earth frame
    u = acc(i,:).';
    u = quat2CEB(quat(i,:)) * u;
    u(3) =  u(3) - 9.81 - x(3);
    
%     % assume the bias affects only acc z in body frame
%     u = acc(i,:).';
%     u(3) = u(3) - x(3);
%     u = quat2CEB(quat(i,:)) * u;
%     u(3) = u(3) - 9.81;
    
    % prediction (this should be done at a constant, fast rate)
%     x = x + ...
%         Ts*[x(2) + 0.5*u(3)*Ts; ...
%             u(3); ...
%             -wa*x(3)] - ...
%         K * x(1);
    x = x + ...
        Ts*[x(2) + 0.5*u(3)*Ts; ...
            u(3); ...
            -wa*x(3)];
     
    % correction (this could also be done slower and not continously, you need to scale the gain K accordingly, lets say every second update: then its 2*K)
%     if mod(downsample_cntr, downsample_divider) == 0
% %         x = x + downsample_divider * K * pos(i);
%         x = x + downsample_divider * K * (pos(i)  - x(1));
%     end
%     downsample_cntr = downsample_cntr + 1;
    if mod(downsample_cntr, downsample_divider) == 0
        downsample_pos = pos(i);
    end
    downsample_cntr = downsample_cntr + 1;
%     x = x + K *downsample_pos;
    x = x + K * (downsample_pos  - x(1));
    
    x_ret(i,:) = x.';
    acc_z_ret(i,1) = u(3);
    acc_z = acc(i,:).';
    acc_z = quat2CEB(quat(i,:)) * acc_z;
    acc_z_ret(i,2) = acc_z(3) - 9.81;

end

end

