clc, clear variables
%%

% notes:
% - samplingrate 50 Hz
% - ind_repair must be set by hand (either 1:2:N or 2:2:N)

is_matlab = true;  % false -> octave
do_plot_raw_signals = false;
do_repair_signals   = false;
do_use_baro = true; % false -> gps

downsample_divider = 1; % if > 1 gps and baro will be downsampled

% 1/Ts/downsample_divider : 1/Ts/downsample_divider : 1/2/Ts

% data = read_bin_data('measurements/20220714/log_40802.bin');
% T_eval = [60 108]; % relative to data.ti
% T_comp = [0 inf];  % relative to T_eval
% data = read_bin_data('measurements/20220714/log_40803.bin');
% T_eval = [1 inf]; % relative to data.ti
% T_comp = [0 inf];  % relative to T_eval
% data = read_bin_data('measurements/20220714/log_40805.bin');
load data_log_40805.mat
T_eval = [1 inf]; % relative to data.ti
T_comp = [0 inf]; % relative to T_eval
baro_offset = 0;

if is_matlab
  addpath ../fcn_bib
else
  pkg load control
  pkg load signal
  pkg load communications
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

Ts = 1/50;
N = size(data.ti, 1);

ind_eval = data.ti >= T_eval(1) & data.ti < T_eval(2);

% extract data so that we are independent from naming above
time = data.ti(ind_eval); time = time - time(1);
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
wa = 2*pi*0.05; % set this zero to get integrator for bias
A = [[0 1 0]; [0 0 -1]; [0 0 -wa]];
B = [0 1 0].';
C = [1 0 0];

% % steady state time continous kalman filter
% Q = diag([1 20 10]);
% R = 0.1*5.0;
% [K, S, e] = lqr(A.', C.', Q, R);
% K = K.';

%%

% prewarp: w0 = 2/Ts*tan(w0*Ts/2);

% pole placement (easy analytical solution)
% % bessel
% w0 = 2*pi*0.2;
% s1 = 0.941600026533207; % 0.149860298638220*2*pi;
% s2 = 1.030544545438434; % 0.164016258482917*2*pi;
% D2 = 0.723540179945206;
% w1 = w0 * s1
% w2 = w0 * s2
% bf pt3
w0 = 2*pi*0.2;
D2 = 1.0;
w1 = w0;
w2 = w0;
% % arbitary (w1, w2, D2)
% w1 = 2*pi*0.1136;
% w2 = 2*pi*0.3983; D2 = 0.7480;
% K = place(A.', C.', [-w1, -w2*D2 + 1i*w2*sqrt(1 - D2^2), -w2*D2 - 1i*w2*sqrt(1 - D2^2)]).'
% analytical solution for place
k1 = (w1 + 2*D2*w2) - wa;
k2 = (w2^2 + 2*D2*w1*w2) - k1*wa;
k3 = k2*wa - w1*w2^2;
K = [k1, k2, k3].';

% % arbitary (w1, w2, w3)
% w1 = 2*pi*0.05;
% w2 = 2*pi*0.1;
% w3 = 2*pi*0.2;
% % K = place(A.', C.', [-w1, -w2, -w3]).'
% % analytical solution for place
% k1 = w1 + w2 + w3 - wa;
% k2 = (w1*w2 + w1*w3 + w2*w3) - k1*wa;
% k3 = k2*wa - w1*w2*w3;
% K = [k1, k2, k3].';

% choose the input and output signal, for the linear filter the bias is on acc_z in earth frame
if do_use_baro
    u = [acc_earth(:,3) - 9.81, baro];
else
    u = [acc_earth(:,3) - 9.81, gps_pos(:,3)];
end

% closed loop system
sys = ss(A - K*C, [B, K], eye(3), 0);
y = lsim(sys, u, time);

% the sum of these give the position and the velocity estimate
sys_inp  = ss(A - K*C, B, C, 0);        % acts on acc
sys_out  = ss(A - K*C, K, C, 0);        % acts on pos
sys_dinp = ss(A - K*C, B, [0 1 0], 0);
sys_dout = ss(A - K*C, K, [0 1 0], 0);

s = tf('s');

% they are the complementary parts (sum up to 1)
G_inp = minreal(tf(sys_inp));
G_out = minreal(tf(sys_out));
G_hp = G_inp * s * s;
G_lp = G_out;

G_dinp = minreal(tf(sys_dinp));
G_dout = minreal(tf(sys_dout));
G_dhp = G_inp * s;
G_dlp = G_out / s;

% % reconstructing those is a pain, these only hold if wa = 0
% tau = 1/(2*pi*0.3);
% G_hp_ = tf([tau^3 0 0 0], conv([tau 1], conv([tau 1], [tau 1])));
% G_lp_ = 1 - G_hp_;
% ---
% G_hp_ = s^3 * tf(1, [1/w1 1])* tf(w2^2, [1 2*D2*w2 w2^2]) /( w1 * w2^2);
% G_lp_ = 1 - G_hp_;
% G_out_ = G_lp_;
% G_inp_ = s * tf(1, [1/w1 1])* tf(w2^2, [1 2*D2*w2 w2^2]) /( w1 * w2^2); % G_hp_ / Gd / Gd;
% ---
% w = w1;
% G_lp_ = tf(1, [1/w 1]) * tf(1, [1/w 1]) * tf(1, [1/w 1]);
% G_hp_ =  tf([1/w 0], [1/w 1]) * tf([1/w 0], [1/w 1]) * tf([1/w 0], [1/w 1]);
% G_lp_ = tf(1, [1/w1 1]) * tf(w2^2, [1 2*D2*w2  w2^2]);
% G_hp_ =  s^3 * tf(1, [1/w1 1])* tf(w2^2, [1 2*D2*w2 w2^2]) /( w1 * w2^2);

% this is exactly the same like y if wa = 0
y_cf = lsim(sys_inp, u(:,1), time) + lsim(sys_out, u(:,2), time);
y_cf = [y_cf, lsim(sys_dinp, u(:,1), time) + lsim(sys_dout, u(:,2), time)];

figure(13)
bode(sys_inp, sys_out, sys_dinp, sys_dout), grid on%, xlim([1e-3 1/2/Ts])
title('time continous filters')
legend('filter for acc \rightarrow pos', 'filter for baro/gps \rightarrow pos', ...
       'filter for acc \rightarrow vel', 'filter for baro/gps \rightarrow vel', 'location', 'northeast')

figure(14)
title('time continous filters')
bode(G_hp, G_lp, G_hp + G_lp), grid on%, xlim([1e-3 1/2/Ts])
legend('highpass', 'lowpass', 'sum of both', 'location', 'northeast')

figure(15)
title('time continous filters')
bode(G_dhp, G_dlp, G_dhp + G_dlp), grid on%, xlim([1e-3 1/2/Ts])
legend('highpass', 'lowpass', 'sum of both', 'location', 'northeast')

%

Gf = c2d(tf(1, [1/(2*pi*4) 1]), Ts, 'tustin');
dpos = diff(u(:,2)) ./ diff(time); dpos = [dpos; dpos(end)];
dlidar = diff(lidar) ./ diff(time); dlidar = [dlidar; dlidar(end)];

% nd_pos corresponds to posDiscreteDelay
nd_pos = 7;
[y_est, acc_z_est] = altitude_estimator(K, nd_pos, wa, acc, u(:,2), quat, Ts);
format long
single(y_est(1:10,:))
format short

% % euler discretization
% Ad = Ts*A + eye(size(A));
% Bd = Ts*B;
% Cd = eye(3);
% Dd = zeros(3,1);

% tustin
Ad = [[1 Ts -Ts^2/2]; ...
      [0 1  -Ts    ];...
      [0 0   1     ]];
Bd = [Ts^2/2, Ts, 0].';
% Cd = [[1 Ts/2 -Ts^2/4]; ...
%       [0 1    -Ts/2  ]; ...
%       [0 0     1]];
% Dd = [Ts^2/4, Ts/2, 0].';
Cd = eye(3);
Dd = zeros(3,1);

% build the discrete time linear filter
sys_inpd  = ss(Ad - K*Ts*Cd(1,:), Bd - K*Ts*Dd(1,1), Cd(1,:), Dd(1,1), Ts);
sys_outd  = ss(Ad - K*Ts*Cd(1,:), K*Ts             , Cd(1,:), Dd(1,1), Ts);
sys_dinpd = ss(Ad - K*Ts*Cd(1,:), Bd - K*Ts*Dd(1,1), Cd(2,:), Dd(2,1), Ts);
sys_doutd = ss(Ad - K*Ts*Cd(1,:), K*Ts             , Cd(2,:), Dd(2,1), Ts);

figure(16)
bode(sys_inpd, sys_outd, sys_dinpd, sys_doutd), grid on%, xlim([1e-3 1/2/Ts])
title('time discrete filters')
legend('filter for acc \rightarrow pos', 'filter for baro/gps \rightarrow pos', ...
       'filter for acc \rightarrow vel', 'filter for baro/gps \rightarrow vel', 'location', 'northeast')

if is_matlab
    % build the discrete time linear filter
    [aa, bb, cc, dd] = dlinmod('altitude_estimator_G', Ts);
    sysd = ss(aa, bb, cc, dd, Ts);

    figure(17)
    bode(sysd(1,1), sysd(1,2), sysd(2,1), sysd(2,2)), grid on%, xlim([1e-3 1/2/Ts])
    title('time discrete filters with delay correction')
    legend('filter for acc \rightarrow pos', 'filter for baro/gps \rightarrow pos', ...
           'filter for acc \rightarrow vel', 'filter for baro/gps \rightarrow vel', 'location', 'northeast')
end

%%

figure(18)
plot(time, u(:,2), 'k'), grid on, hold on
plot(time, y(:,1), 'Linewidth', 2, 'color', [0 0 1])
plot(time, y_est(:,1), 'Linewidth', 2, 'color', [0 0.5 0])
plot(time, y_cf(:,1), 'c', 'Linewidth', 2)
plot(time, lidar, 'Linewidth', 2, 'color', [1 0 0])
plot(time, est_xyz(:,3), 'm', 'Linewidth', 2)
plot(time, est_z_2, 'Linewidth', 2, 'color', [0.6 0.6 0.6]), hold off
ylabel('Pos z (m)'), xlabel('Time (s)'), xlim([0 time(end)])%, ylim([-1 7])

figure(19)
plot(time, filtfilt(Gf.num{1}, Gf.den{1}, dpos), 'k'), grid on, hold on
plot(time, y(:,2), 'Linewidth', 2, 'color', [0 0 1])
plot(time, y_est(:,2), 'Linewidth', 2, 'color', [0 0.5 0])
plot(time, y_cf(:,2), 'c', 'Linewidth', 2)
plot(time, filtfilt(Gf.num{1}, Gf.den{1}, dlidar), 'Linewidth', 2, 'color', [1 0 0])
plot(time, est_Vxyz(:,3), 'm', 'Linewidth', 2)
plot(time, est_Vz_2, 'Linewidth', 2, 'color', [0.6 0.6 0.6]), hold off
ylabel('Vel z (m/s)'), xlabel('Time (s)'), xlim([0 time(end)])%, ylim([-5 5])

figure(20)
plot(time, [y(:,3), y_est(:,3)], 'Linewidth', 2), grid on, hold off
ylabel('Acc Bias (m/s^2)'), xlabel('Time (s)'), xlim([0 time(end)])

figure(21)
plot(time, acc_z_est, 'Linewidth', 2), grid on, xlim([0 time(end)])
ylabel('Acc without Bias (m/s^2)'), xlabel('Time (s)'), xlim([0 time(end)])

%%

ind_comp = T_comp(1) <= time & time <= T_comp(2);

finddelay(y_est(ind_comp,1), est_xyz(ind_comp,3))

Nest     = round(5/Ts);
koverlap = 0.95;
Noverlap = round(koverlap*Nest);
window   = hann(Nest);
delta    = 0;
[G1, C1, freq] = my_tfestimate(acc_z_est(ind_comp,1), y_est(ind_comp,1), window, Noverlap, Nest, Ts, delta);
[G2, C2, freq] = my_tfestimate(y_est(ind_comp,1), est_xyz(ind_comp,3), window, Noverlap, Nest, Ts, delta);

figure(21)
if is_matlab
  bode(G1, G2), grid on, xlim([G1.Frequency(2) 1/2/Ts])
else
  ind = G1.w > 0;
  G = [squeeze(G1.H), squeeze(G2.H)];
  subplot(211)
  loglog(G1.w(ind), abs(G(ind,:))), grid on
  ylabel('Magnitude (abs)')
  subplot(212)
  semilogx(G1.w(ind), 180/pi*angle(G(ind,:))), grid on
  ylabel('Phase (deg)')
end

% Nest     = size(u, 1); %round(2.5/Ts);
% koverlap = 0;
% Noverlap = round(koverlap*Nest);
% window   = ones(Nest, 1); % hann(Nest);
% data_spectras = y_est(:,1);
% [pxx, f] = pwelch(data_spectras, window, Noverlap, 1/Ts, 'Power');
% [pxx, f] = pwelch(data_spectras, window, Noverlap, Nest, 1/Ts, 'Power');
% spectras = sqrt(pxx*2); % power -> amplitude (dc needs to be scaled differently)
[P_est, f]   = my_fft(y_est(:,1), Ts);
[P_cf, f]    = my_fft(y_cf(:,1), Ts);
[Vel_est, f] = my_fft(y_est(:,2), Ts);
[Acc_est, f] = my_fft(acc_z_est, Ts);
ind = f > 0;

figure(22)
subplot(311)
plot(f(ind), [P_est(ind), P_cf(ind)]), grid on, xlabel('Frequency (Hz)')
set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log')
subplot(312)
plot(f(ind), Vel_est(ind)), grid on, xlabel('Frequency (Hz)')
set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log')
subplot(313)
plot(f(ind), Acc_est(ind)), grid on, xlabel('Frequency (Hz)')
set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log')

%% convert pi-d to p-pi

% clc, clear all
% syms kp ki kd Kv P I
% eqns = [kp == Kv*P + I, ki == Kv*I, kd == P];
% vars = [Kv, P, I];
% sol = solve(eqns, vars)
% % sol.Kv
% (kp/2 + (kp^2 - 4*kd*ki)^(1/2)/2)/kd
% (kp/2 - (kp^2 - 4*kd*ki)^(1/2)/2)/kd
% % sol.P
% kd
% kd
% % sol.I
% kp/2 - (kp^2 - 4*kd*ki)^(1/2)/2
% kp/2 + (kp^2 - 4*kd*ki)^(1/2)/2

%% pole placemant K

% clc, clear all
% syms  wa k1 k2 k3 s w1 w2 D2 w3
% A = [[0 1 0]; [0 0 -1]; [0 0 -wa]]
% C = [1 0 0]
% K = [k1; k2; k3]
% collect(det(s*eye(3) - (A - K*C)), 's')
% % -> s^3 + (k1 + wa)*s^2 + (k2 + k1*wa)*s - k3 + k2*wa
%
% collect(w1*(1/w1*s + 1)*(s^2 + 2*D2*w2*s + w2^2), 's')
% % -> s^3 + (w1 + 2*D2*w2)*s^2 + (w2^2 + 2*D2*w1*w2)*s + w1*w2^2
%
% collect((w1*w2*w3)*(1/w1*s + 1)*(1/w2*s + 1)*(1/w3*s + 1), 's')
% % -> s^3 + (w1 + w2 + w3)*s^2 + (w1*w2 + w1*w3 + w2*w3)*s + w1*w2*w3

%%

% clc, clear all
% syms w1 w2 w3 s
% f = (w1 + w2 + w3)*s^2 + (w1*w2 + w1*w3 + w2*w3)*s + w1*w2*w3
% sol = simplify(solve(f, s))
