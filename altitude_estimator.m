function [x_ret, acc_z_ret] = altitude_estimator(K, wa, acc, pos, quat, Ts)

N = size(acc, 1);
x = zeros(3,1);
x_ret = zeros(N,3);
acc_z_ret = zeros(N,1);

% p = p + v*Ts + 0.5*a*Ts^2;
% v = v + a*Ts;
% b = b;
for i = 1:N
    
%     % directly take acc
%     u = acc(i,:).';
%     u(3) = u(3) - 9.81 - x(3);

%     % assume the bias affects only acc z in earth frame
%     u = acc(i,:).';
%     u = quat2CEB(quat(i,:)) * u;
%     u(3) =  u(3) - 9.81 - x(3);
    
    % assume the bias affects only acc z in body frame
    u = acc(i,:).';
    u(3) = u(3) - x(3);
    u = quat2CEB(quat(i,:)) * u;
    u(3) = u(3) - 9.81;
    
    % prediction (this should be done at a constant, fast rate)
    x = x + ...
        Ts*[x(2) + 0.5*u(3)*Ts; ...
            u(3); ...
            -wa*x(3)] - ...
        K * x(1);
     
    % correction (this could also be done slower and not continously, you need to scale the gain K accordingly, lets say every second update: then its 2*K)
    x = x + K * pos(i);
    
    x_ret(i,:) = x.';
    acc_z_ret(i) = u(3);

end

end

