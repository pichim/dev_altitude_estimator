clc, clear variables
%%

% notes:
% - samplingrate 50 Hz
% - ind_repair must be set by hand (either 1:2:N or 2:2:N)

is_matlab = true; % false -> octave
do_plot_raw_signals = false;

load data_log_40802.mat

if is_matlab
  addpath ../fcn_bib
else
  pkg load control
  pkg load signal
  addpath fcn_bib_copy
end
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

T_eval = [8 261];
Ts = 1/50;
N = size(data.ti, 1);
baro_offset = 443.9;

ind_eval = data.ti >= T_eval(1) & data.ti < T_eval(2);

% extract data so that we are independent from naming above
time = data.ti(ind_eval); time = time - time(1);
acc = data.acc(ind_eval,:);
est_rpy = data.est_RPY(ind_eval,:);
est_xyz = data.est_xyz(ind_eval,:);

ind_repair = data.Lidar(:,4) > 0; % repair lidar and append as last coloumn
lidar = interp1(data.ti(ind_repair), data.Lidar(ind_repair,3), data.ti, 'linear', 'extrap');
lidar = lidar(ind_eval,:);

% gps_pos = data.GPS_pos(ind_eval,:);
% gps_vel = data.GPS_vel(ind_eval,:);
%ind_first = find(data.GPS_pos(:,1) ~= 0, 1); % this is just a helper
ind_repair = 1:2:N;
gps_pos = interp1(data.ti(ind_repair), data.GPS_pos(ind_repair,:), data.ti, 'linear', 'extrap');
gps_pos = gps_pos(ind_eval,:);
gps_vel = interp1(data.ti(ind_repair), data.GPS_vel(ind_repair,:), data.ti, 'linear', 'extrap');
gps_vel = gps_vel(ind_eval,:);

ind_repair = 1:2:N;
baro = interp1(data.ti(ind_repair), data.Baro(ind_repair,1), data.ti, 'linear', 'extrap');
baro = baro(ind_eval,1) - baro_offset;

quat = rpy2quat(data.est_RPY);
acc_earth = zeros(N,3);
for i = 1:N
    CEB = quat2CEB(quat(i,:));
    acc_earth(i,:) = ( CEB * data.acc(i,:).' ).';
end
quat = quat(ind_eval,:);
acc_earth = acc_earth(ind_eval,:);

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

% time continous model
wa = 0*2*pi*0.02; % set this zero to get integrator for bias
A = [[0 1 0]; [0 0 -1]; [0 0 -wa]];
B = [0 1 0].';
C = [1 0 0];

% euler forward discretization, you need to use dlqr for a discrete time system
% A = Ts*A + eye(size(A));
% B = Ts*B;

% static kalman filter
% Q = diag([1 1 1e0]);
% R = 1e1;
% [K, S, e] = lqr(A.', C.', Q, R);

% pole placement (easy analytical solution)
% w1 <= w2, 0.6 <= D <= 1 (D = 1/2/Q)
w1 = 2*pi*0.15;
w2 = w1; D2 = 0.8;
% K = place(A.', C.', [-w1, -w2*D2 + 1i*w2*sqrt(1 - D2^2), -w2*D2 - 1i*w2*sqrt(1 - D2^2)]).'
% analytical solution from place
k1 = (w1 + 2*D2*w2) - wa;
k2 = (w2^2 + 2*D2*w1*w2) - k1*wa;
k3 = k2*wa - w1*w2^2;
K = [k1, k2, k3].';

% choose the input and output signal
u = [acc(:,3) - 9.81, baro]; % acc_earth(:,3)

% closed loop system
sys = ss(A - K*C, [B, K], eye(3), 0);
y = lsim(sys, u, time);

% since we calculated K based on time continous model we need to scale it with Ts
[y_est, acc_z_est] = altitude_estimator(K*Ts, wa, acc, u(:,2), quat, Ts);

% the sum of these give the position and the velocity estimate
sys_inp = ss(A - K*C, B, C, 0);        % acts on acc
sys_out = ss(A - K*C, K, C, 0);        % acts on pos
sys_dinp = ss(A - K*C, B, [0 1 0], 0);
sys_dout = ss(A - K*C, K, [0 1 0], 0);

s = tf('s');

% % the are the complementary parts (sum up to 1)
% G_inp = minreal(tf(sys_inp));
% G_out = minreal(tf(sys_out));
% G_hp = G_inp * s * s;
% G_lp = G_out;

% reconstructing those is a pain, these only hold if wa = 0
G_hp = s^3 * tf(1, [1/w1 1])* tf(w2^2, [1 2*D2*w2 w2^2]) /( w1 * w2^2);
G_lp = 1 - G_hp;
G_out = G_lp;
G_inp = s * tf(1, [1/w1 1])* tf(w2^2, [1 2*D2*w2 w2^2]) /( w1 * w2^2); % G_hp / Gd / Gd;

% this is exactly the same like y if wa = 0
y_cf = lsim(G_inp, u(:,1), time) + lsim(G_out, u(:,2), time);

figure(13)
bode(G_inp, G_out), grid on%, xlim([1e-3 1/2/Ts])
legend('filter for acc', 'filter for pos', 'location', 'best')

figure(14)
bode(G_hp, G_lp, G_hp + G_lp), grid on%, xlim([1e-3 1/2/Ts])
legend('highpass', 'lowpass', 'sum of both', 'location', 'best')

Gf = c2d(tf(1, [1/(2*pi*2) 1]), Ts, 'tustin');
dpos = diff(u(:,2)) ./ diff(time); dpos = [dpos; dpos(end)];
dlidar = diff(lidar) ./ diff(time); dlidar = [dlidar; dlidar(end)];

figure(15)
plot(time, u(:,2), 'k'), grid on, hold on
plot(time, y(:,1), 'Linewidth', 2, 'color', [0 0 1])
plot(time, y_est(:,1), 'Linewidth', 2, 'color', [0 0.5 0])
plot(time, lidar, 'Linewidth', 2, 'color', [1 0 0]), hold off
ylabel('Pos z (m)'), xlabel('Time (s)'), xlim([0 time(end)])

figure(16)
plot(time, filtfilt(Gf.num{1}, Gf.den{1}, dpos), 'k'), grid on, hold on
plot(time, y(:,2), 'Linewidth', 2, 'color', [0 0 1])
plot(time, y_est(:,2), 'Linewidth', 2, 'color', [0 0.5 0])
plot(time, filtfilt(Gf.num{1}, Gf.den{1}, dlidar), 'Linewidth', 2, 'color', [1 0 0]), hold off
ylabel('Vel z (m/s)'), xlabel('Time (s)'), xlim([0 time(end)])

figure(17)
plot(time, [y(:,3), y_est(:,3)], 'Linewidth', 2), grid on, hold off
ylabel('Acc Bias (m/s^2)'), xlabel('Time (s)'), xlim([0 time(end)])

Nest     = round(2/Ts);
koverlap = 0.95;
Noverlap = round(koverlap*Nest);
window   = hann(Nest);
delta    = 0;
[G1, C1, freq] = my_tfestimate(acc_z_est, y_est(:,1), window, Noverlap, Nest, Ts, delta);
[G2, C2, freq] = my_tfestimate(y_est(:,1), lidar, window, Noverlap, Nest, Ts, delta);

figure(18)
if is_matlab
  bode(G1, G2), grid on, xlim([G1.Frequency(2) 1/2/Ts])
else
  subplot(211)
  loglog(2*pi*G1.w, abs(squeeze([G1.H, G2.H]))), grid on
  ylabel('Magnitude (abs)')
  subplot(212)
  semilogx(2*pi*G1.w, 180/pi*angle(squeeze([G1.H, G2.H]))), grid on
  ylabel('Phase (deg)')
end

%% some symbolic stuff to calculate the pole placemant K

% clc, clear all
% syms  wa k1 k2 k3 s w1 w2 D2
% A = [[0 1 0]; [0 0 -1]; [0 0 -wa]]
% C = [1 0 0]
% K = [k1; k2; k3]
% collect(det(s*eye(3) - (A - K*C)), 's')
% % -> s^3 + (k1 + wa)*s^2 + (k2 + k1*wa)*s - k3 + k2*wa
%
% collect(w1*(1/w1*s + 1)*(s^2 + 2*D2*w2*s + w2^2), 's')
% % -> s^3 + (w1 + 2*D2*w2)*s^2 + (w2^2 + 2*D2*w1*w2)*s + w1*w2^2
