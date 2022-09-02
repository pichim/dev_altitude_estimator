function [x_ret, acc_z_ret] = position_estimator(f_cut, f_a, Ts, nd_pos, acc, pos, quat)

N = size(acc, 1);
position = 0;
velocity = 0;
accBias  = 0;
x_ret = zeros(N,3);
acc_z_ret = zeros(N,1);

pos_past_idx = 1;
pos_past = zeros(1+nd_pos, 1); % ringbuffer

k  = Ts/(1/(2*pi*f_cut) + Ts);
wa = 2*pi*f_a;
a33 = 1/(Ts*wa + 1);

% Ad = [[1, Ts,  -Ts^2*a33]
%       [0,  1,  -Ts*a33]
%       [0,  0,  a33]];   
% Bd = [Ts^2
%       Ts
%       0];

% define location of discrete poles
w1 = (k-1);
w2 = w1;
w3 = w1;
% calculate observer gain
c1 = 1 - (w1*w2 + w1*w3 + w2*w3);
k1 = (w1 * w2 * w3) / a33 + 1;
c2 = (c1 - k1) / (Ts*a33);
k2 = c2 + (2 - k1)/Ts;
k3 = (c2/a33 - ((w1 + w2 + w3)/a33 + 1)/Ts)/Ts;

for i = 1:N
    
    % CEB =
    % [cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)]
    % [cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi)]
    % [        -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta)]
    
    
    % assume the bias affects only acc z in body frame
    acc_i = acc(i,:).';
    acc_i(3) = acc_i(3) - a33 * accBias;
    acc_i = quat2CEB(quat(i,:)) * acc_i;
    acc_i = acc_i(3) - 9.81;
    
%     % assume the bias affects only acc z in earth frame
%     acc_i = acc(i,:).';
%     acc_i = quat2CEB(quat(i,:)) * acc_i;
%     acc_i = acc_i(3) - 9.81 - a33 * accBias;
        
    % update delayed position estimation
    % pos_past(pos_past_idx) = Ad(1,:) * x + Bd(1,:) * acc_i;
    pos_past(pos_past_idx) = position + Ts*(velocity + Ts*acc_i);
    pos_past_idx = pos_past_idx + 1;
    if pos_past_idx > 1+nd_pos
        pos_past_idx = 1;
    end
    
    %  x = Ad * x + Bd * acc_i + K * (pos(i) - pos_past(pos_past_idx));
    position = position + Ts*(velocity + Ts*acc_i);
    velocity = velocity + Ts*acc_i ;
    accBias = a33 * accBias;
    
    estimationError = pos(i) - pos_past(pos_past_idx);
    position = position + k1 * estimationError;
    velocity = velocity + k2 * estimationError;
    accBias  = accBias  + k3 * estimationError;
           
    x_ret(i,:) = [position, velocity, accBias];
    acc_z_ret(i,1) = acc_i;

end

end

