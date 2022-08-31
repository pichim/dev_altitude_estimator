clc, clear all
%%

try
  addpath ../fcn_bib
catch
  addpath fcn_bib_copy
end

% notes:
% - samplingrate 50 Hz
% - ind_repair must be set by hand (either 1:2:N or 2:2:N)

do_plot_raw_signals = false;
do_repair_signals   = false;
do_use_baro = true; % false -> gps

downsample_divider = 1; % if > 1 gps and baro will be downsampled

% 1/Ts/downsample_divider : 1/Ts/downsample_divider : 1/2/Ts


data = read_bin_data('measurements/20220708/log_40807.bin');
T_eval = [1 inf]; % relative to data.ti
T_comp = [0 inf]; % relative to T_eval
baro_offset = -0.55;

% data = read_bin_data('measurements/20220711/log_40803.bin');
% T_eval = [2.02+13.54 192]; % relative to data.ti
% T_comp = [0 inf]; % relative to T_eval
% baro_offset = 376.93;

% data = read_bin_data('measurements/20220712/log_40804.bin');
% T_eval = [1 inf]; % relative to data.ti
% T_comp = [0 inf]; % relative to T_eval
% baro_offset = 0;

% data = read_bin_data('measurements/20220714/log_40802.bin');
% T_eval = [60 108]; % relative to data.ti
% T_comp = [0 inf]; % relative to T_eval
% baro_offset = 0;

% % data = read_bin_data('measurements/20220714/log_40803.bin');
% % T_eval = [1 inf]; % relative to data.ti
% % T_comp = [0 inf]; % relative to T_eval
% % baro_offset = 0;

% data = read_bin_data('measurements/20220714/log_40805.bin');
% T_eval = [1 inf]; % relative to data.ti
% T_comp = [0 inf]; % relative to T_eval
% baro_offset = 0;

% data =
%   struct with fields:
%                ti: [13085×1  double]
%           alldata: [13085×63 double]
%               gyr: [13085×3  double]
%               acc: [13085×3  double]
%                RS: [13085×7  double]
%            RS_yaw: [13085×1  double]
%           est_RPY: [13085×3  double]
%           est_xyz: [13085×3  double]
%        cntrl_Mxyz: [13085×3  double]
%          cntrl_FT: [13085×1  double]
%     cntrl_xyz_des: [13085×3  double]
%             sm_FM: [13085×1  double]
%             Lidar: [13085×5  double]
%             sm_FS: [13085×1  double]
%     cntrl_RPY_des: [13085×3  double]
%          est_Vxyz: [13085×3  double]
%           GPS_pos: [13085×3  double] 25 Hz
%           GPS_vel: [13085×3  double] 25 Hz
%          GPS_head: [13085×2  double] 25 Hz
%        GPS_numSat: [13085×1  double] 25 Hz
%          GPS_time: [13085×1  double] 25 Hz
%                OF: [13085×4  double] 25 Hz
%               mag: [13085×3  double]
%              mag2: [13085×3  double]
%              Baro: [13085×3  double] 25 Hz

%%

Ts = 1/50;
N = size(data.ti, 1);

ind_eval = data.ti >= T_eval(1) & data.ti < T_eval(2);

% extract data so that we are independent from naming above
time = data.ti(ind_eval); time = time - time(1);

% data.acc(:,3) = data.acc(:,3) + 5;

acc = data.acc(ind_eval,:);
est_rpy = data.est_RPY(ind_eval,:);
est_xyz = data.est_xyz(ind_eval,:);
est_Vxyz = data.est_Vxyz(ind_eval,:);
est_z_2 = data.est_z_2(ind_eval,:) - baro_offset;
est_Vz_2 = data.est_Vz_2(ind_eval,:);
cntrl_FT = data.cntrl_FT(ind_eval,:);

if do_repair_signals
    ind_repair = data.Lidar(:,4) > 0; % repair lidar and append as last coloumn
    lidar = interp1(data.ti(ind_repair), data.Lidar(ind_repair,3), data.ti, 'linear', 'extrap');
    lidar = lidar(ind_eval,:);
else
    lidar = data.Lidar(ind_eval,3);
end

if do_repair_signals
    ind_repair = 1:2:N;
    gps_pos = interp1(data.ti(ind_repair), data.GPS_pos(ind_repair,:), data.ti, 'linear', 'extrap');
    gps_pos = gps_pos(ind_eval,:);
    gps_vel = interp1(data.ti(ind_repair), data.GPS_vel(ind_repair,:), data.ti, 'linear', 'extrap');
    gps_vel = gps_vel(ind_eval,:);
else
    gps_pos = data.GPS_pos(ind_eval,:);
    gps_vel = data.GPS_vel(ind_eval,:);
end

if do_repair_signals
    ind_repair = 1:2:N;
    baro = interp1(data.ti(ind_repair), data.Baro(ind_repair,1), data.ti, 'linear', 'extrap');
    baro = baro(ind_eval,1) - baro_offset;
else
    baro = data.Baro(ind_eval,1) - baro_offset;
end

if downsample_divider > 1
    downsample_cntr = 0;
    for i = 1:size(gps_pos, 1)
        if mod(downsample_cntr, downsample_divider) == 0
            gps_pos_i = gps_pos(i,:);
            gps_vel_i = gps_vel(i,:);
            baro_i = baro(i,:);
        end
        gps_pos(i,:) = gps_pos_i;
        gps_vel(i,:) = gps_vel_i;
        baro(i,:) = baro_i;
        downsample_cntr = downsample_cntr + 1;
    end
end

quat = rpy2quat(data.est_RPY);
acc_earth = zeros(N,3);
for i = 1:N
    CEB = quat2CEB(quat(i,:));
    acc_earth(i,:) = ( CEB * data.acc(i,:).' ).';
end
quat = quat(ind_eval,:);
acc_earth = acc_earth(ind_eval,:);

N = size(time, 1);

%%

if do_plot_raw_signals

    figure(1)
    plot(diff(data.ti(ind_eval))), grid on

    figure(2)
    plot(time, data.gyr(ind_eval,:)), grid on
    ylabel('gyr (rad/s)'), xlabel('Time (sec)')

    figure(3)
    plot(time, [acc, acc_earth]), grid on
    ylabel('acc (m/s^2)'), xlabel('Time (sec)')

    figure(4)
    plot(time, data.RS(ind_eval,:)), grid on
    ylabel('RS ()'), xlabel('Time (sec)')

    figure(5)
    plot(time, data.RS_yaw(ind_eval,:) * 180/pi), grid on
    ylabel('RS yaw (deg)'), xlabel('Time (sec)')

    figure(6)
    plot(time, est_rpy * 180/pi), grid on
    ylabel('est RPY (deg)'), xlabel('Time (sec)')

    figure(7)
    plot(time, est_xyz), grid on
    ylabel('est xyz (m)'), xlabel('Time (sec)')

    figure(8)
    subplot(211)
    plot(time, data.Lidar(ind_eval,[1 3])), grid on
    ylabel('Lidar (m)'), xlabel('Time (sec)')
    subplot(212)
    plot(time, lidar), grid on
    ylabel('Lidar (m)'), xlabel('Time (sec)')

    figure(9)
    plot(time, gps_pos), grid on
    ylabel('gps pos (m)'), xlabel('Time (sec)')

    Gf = tf([1 -1], [Ts 0], Ts);

    figure(10)
    subplot(211)
    plot(time, gps_vel), grid on
    ylabel('gps vel (m/s)'), xlabel('Time (sec)')
    subplot(212)
    plot(time, [gps_vel, filter(Gf.num{1}, Gf.den{1}, gps_pos)]), grid on
    ylabel('gps vel f (m/s)'), xlabel('Time (sec)')

    figure(11)
    plot(time, baro), grid on
    ylabel('baro (m)'), xlabel('Time (sec)')

    figure(12)
    plot(time, [baro, lidar, gps_pos(:,3)]), grid on

end

%%

f_cut = 0.1;
wa = 2*pi*0.01; % this will be beneficial for horizontal position estimates (x-y)
nd_pos = 3;      % this is somewhat a hack and can lead to a unstable filter

% choose the input and output signal, for the linear filter the bias is on acc_z in earth frame
if do_use_baro
    u = [acc_earth(:,3) - 9.81, baro];
else
    u = [acc_earth(:,3) - 9.81, gps_pos(:,3)];
end

%%

k  = Ts/(1/(2*pi*f_cut) + Ts);

a33 = 1/(Ts*wa + 1);

Ad = [[1, Ts,  -Ts^2*a33]
      [0,  1,  -Ts*a33]
      [0,  0,  a33]];
Bd = [Ts^2
      Ts
      0];

% define location of discrete poles
w1 = (k-1);
w2 = w1;
w3 = w1;
% calculate observer gain
c0 = w1*w2*w3;
c1 = 1 - (w1*w2 + w1*w3 + w2*w3);
k1 = c0/a33 + 1;
k2 = (c1 - k1 + (2 - k1)*a33)/(Ts*a33);
k3 = (c1 - k1 - ((w1 + w2 + w3) + a33)*a33)/(Ts^2*a33^2);
K = [k1; k2; k3];

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

% apply some betaflight filtering here
[~, B, A] = get_filter('pt2', 15.0, Ts/3);
acc_f = zeros(3*N, 3);
u_f   = zeros(3*N, 2);
for i = 1:N
    acc_f((i-1)*3+1:i*3,:) = acc(i,:) .* ones(3,3);
    u_f  ((i-1)*3+1:i*3,:) = u  (i,:) .* ones(3,2);
end
acc_f = filter(B, A, acc_f);
u_f   = filter(B, A, u_f  );
acc_f = acc_f(1:3:end,:);
u_f   = u_f  (1:3:end,:);
         
% this is the same like above only if nd_pos = 0
[aa, bb, cc, dd] = dlinmod('position_estimator_G');
sys = ss(aa, bb, cc, dd, Ts);

% y_est = dlsim(aa, bb, cc, dd, u);
x0 = [0; 0; 0*9.81];
y_est = lsim(sys, u_f, time, [x0; zeros(nd_pos, 1)]);


[y_est_nl, acc_z_est] = position_estimator(Ts, G, nd_pos, x0, acc_f, u_f(:,2), quat);
format long
G
single(y_est_nl(1:10,:))
format short

%%

figure(1)
plot(time, u(:,2), 'k'), grid on, hold on
plot(time, lidar, 'r', 'Linewidth', 2)
plot(time, est_xyz(:,3), 'c', 'Linewidth', 2)
plot(time, y_est(:,1), 'b', 'Linewidth', 2)
plot(time, y_est_nl(:,1), 'Linewidth', 2, 'color', [0 0.5 0]), hold off
ylabel('Pos z (m)'), xlabel('Time (s)'), xlim([0 time(end)])

dpos = diff(u(:,2)) ./ diff(time); dpos = [dpos; dpos(end)];
dlidar = diff(lidar) ./ diff(time); dlidar = [dlidar; dlidar(end)];
% Gf = c2d(tf(1, [1/(2*pi*4) 1]), Ts, 'tustin');
% dpos   = filtfilt(Gf.num{1}, Gf.den{1}, dpos);
% dlidar = filtfilt(Gf.num{1}, Gf.den{1}, dlidar);

figure(2)
plot(time, dpos, 'k'), grid on, hold on
plot(time, dlidar, 'r', 'Linewidth', 2)
plot(time, est_Vxyz(:,3), 'c', 'Linewidth', 2)
plot(time, y_est(:,2), 'b', 'Linewidth', 2)
plot(time, y_est_nl(:,2), 'Linewidth', 2, 'color', [0 0.5 0]), hold off
ylabel('Vel z (m/s)'), xlabel('Time (s)'), xlim([0 time(end)])%, ylim([-5 5])

figure(3)
plot(time, [y_est(:,3), y_est_nl(:,3)], 'Linewidth', 2), grid on, hold off
ylabel('Acc Bias (m/s^2)'), xlabel('Time (s)'), xlim([0 time(end)])

figure(4)
plot(time, acc_z_est, 'Linewidth', 2), grid on, xlim([0 time(end)])
ylabel('Acc without Bias (m/s^2)'), xlabel('Time (s)'), xlim([0 time(end)])

%%

ind_comp = T_comp(1) <= time & time <= T_comp(2);

finddelay(y_est_nl(ind_comp,1), est_xyz (ind_comp,3))
finddelay(y_est_nl(ind_comp,2), est_Vxyz(ind_comp,3))

Nest     = round(5/Ts);
koverlap = 0.95;
Noverlap = round(koverlap*Nest);
window   = hann(Nest);
delta    = 0;
[G1, C1, freq] = my_tfestimate(acc_z_est(ind_comp,1), y_est_nl(ind_comp,1), window, Noverlap, Nest, Ts, delta);
[G2, C2, freq] = my_tfestimate(y_est_nl(ind_comp,1), est_xyz(ind_comp,3), window, Noverlap, Nest, Ts, delta);

figure(5)
bode(G1, G2, 2*pi*G1.Frequency), grid on, xlim([G1.Frequency(2) 1/2/Ts])
