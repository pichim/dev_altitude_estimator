clc, clear all
%%

try
  addpath ../fcn_bib
catch
  addpath fcn_bib_copy
end

% load data_20220901_00 % save data_20220901_00 data time
% f_cut = 0.1;
% f_a   = 0.0;
% nd_pos = 9;     % this is somewhat a hack and can lead to a unstable filter

load data_20220901_01 % save data_20220901_00 data time
f_cut = 0.2;
f_a   = 0.0;
nd_pos = 9;     % this is somewhat a hack and can lead to a unstable filter

% notes:
% - position.c is running 117.6471 or 125 Hz, not consistent at 120 Hz
% - position.c should run at 125 Hz???

% ind_pid_p   =  3: 5;
% ind_pid_i   =  6: 8;
% ind_pid_d   =  9:10;
% ind_pid_f   = 11:13;
% ind_rc      = 14:16;
% ind_setp    = 18:21;
% ind_baro    = 24;
% ind_gyro_f  = 26:28;
% ind_acc     = 29:31;
% ind_mot     = 36:39;
ind_debug   = 32:35;
% ind_pid_sum = 48:50;
% ind_pid_err = 55:57;
% ind_pi_sum  = 62:64;

% DEBUG_SET(DEBUG_POSITION_ESTIMATOR_Z, 0, accZ);
% DEBUG_SET(DEBUG_POSITION_ESTIMATOR_Z, 1, baroAlt);
% DEBUG_SET(DEBUG_POSITION_ESTIMATOR_Z, 2, positionEstimatorZ.position - baroAltOffset);
% DEBUG_SET(DEBUG_POSITION_ESTIMATOR_Z, 3, positionEstimatorZ.velocity);
% 0: accZ rotated to earth frame
% 1: baro alt filtered and upsampled
% 2: position estimator position
% 3: position estimator velocity

figure(1)
subplot(211)
stairs(time(1:end-1), diff(time))
subplot(212)
stairs(time, data(:, ind_debug(2)))

mean(diff(time))
median(diff(time))

%%

Ts_log = 1/2e3;
Ts_position = 1/125;
ind_eval = 20:Ts_position/Ts_log:size(data, 1);
accZ     = data(ind_eval, ind_debug(1));
baroAlt  = data(ind_eval, ind_debug(2));
position = data(ind_eval, ind_debug(3));
velocity = data(ind_eval, ind_debug(4));
time     = time(ind_eval);

Ts_position = 1/120;

% Ts_log = 1/2e3;
% Ts_position = 1/120;
% accZ_past = inf;
% accZ = [];
% baroAlt = [];
% position = [];
% velocity = [];
% time_copy = time;
% time = [];
% for i = 1:size(data, 1)
%     if data(i, ind_debug(1)) ~= accZ_past
%         accZ_past = data(i, ind_debug(1));
%         accZ(end+1,1)     = data(i, ind_debug(1));
%         baroAlt(end+1,1)  = data(i, ind_debug(2));
%         position(end+1,1) = data(i, ind_debug(3));
%         velocity(end+1,1) = data(i, ind_debug(4));
%         time(end+1,1)     = time_copy(i);
%     end
% end

figure(2)
subplot(221)
stairs(time, accZ), grid on
subplot(222)
stairs(time, baroAlt), grid on
subplot(223)
stairs(time, position), grid on
subplot(224)
stairs(time, velocity), grid on

%%

y_est = position_estimator_debug_bf(f_cut, f_a, nd_pos, Ts_position, accZ, baroAlt);

figure(3)
subplot(211)
stairs(time, [baroAlt, position, y_est(:,1)]), grid on
subplot(212)
stairs(time, velocity, 'color', [0 0.5 0]), grid on, hold on
stairs(time, y_est(:,2), 'r'), hold off

% %%
% 
% N = size(alt_raw, 1);
% 
% alt_f0 = alt_raw;
% n = 9;
% for i = n:N
%     alt_window = alt_raw(i-n+1:i);
%     [~, idx] = sort(alt_window);
%     if mod(n, 2) == 0
%         alt_f0(i) = mean(alt_window(idx(n/2:n/2+1)));
%     else
%         alt_f0(i) = alt_window(idx((n+1)/2));
%     end
% end
% % alt_f0 = filter([0.5 0.5], [1 0], alt_f0);
% 
% alt_f1 = alt_raw;
% % n = 7;
% % for i = n:N
% %     alt_window = alt_raw(i-n+1:i);
% %     [~, idx] = sort(alt_window);
% %     if mod(n, 2) == 0
% %         alt_f1(i) = mean(alt_window(idx(n/2:n/2+1)));
% %     else
% %         alt_f1(i) = alt_window(idx((n+1)/2));
% %     end
% % end
% 
% f_cut = 2.0;
% [G, B, A] = get_filter('pt2', f_cut, Ts_baro);
% alt_f1 = filter(B, A, alt_f1);
% 
% G_fast = get_filter('pt2', f_cut, 1/120);
% figure(99)
% bode(G, G_fast), grid on
% 
% linewidth = 1.0;
% 
% % XLim = [63 68]; YLim = [9200 10600];
% % XLim = [153 161]; YLim = [13200 14600];
% % XLim = [395 400]; YLim = [1400 2600];
% XLim = [25 38]; YLim = [3400 6600];
% % XLim = [0 inf]; YLim = [-inf inf];
% 
% figure(2)
% ax(1) = subplot(221);
% stairs(time, [alt_raw, alt, alt_old], 'Linewidth', linewidth), grid on, ylim(YLim)
% ylabel('Altitude (cm)')
% legend('baroAltitudeRaw (Debug 2)', 'filteredBaroAltitudeRaw (Debug 1)', 'baro.BaroAlt (Debug 3)', 'location', 'best')
% ax(2) = subplot(223);
% stairs(time, throttle, 'k', 'Linewidth', linewidth), grid on
% xlabel('Time (sec)'), ylabel('Throttle')
% ax(3) = subplot(122);
% stairs(time, alt_raw, 'Linewidth', linewidth), grid on, hold on
% stairs(time, alt_f0, 'k', 'Linewidth', linewidth)
% stairs(time, alt_f1, 'm', 'Linewidth', linewidth), hold off, ylim(YLim)
% xlabel('Time (sec)'), ylabel('Altitude (cm)')
% legend('baroAltitudeRaw (Debug 2)', ...
%        'baroAltitudeRaw -> 9 point median (offline)', ...
%        'baroAltitudeRaw -> pt2 filter (offline)', 'location', 'best')
% linkaxes(ax, 'x')
% xlim(XLim)
% 
% figure(3)
% plot(time, alt_raw, 'Linewidth', linewidth), grid on, hold on
% plot(time, alt_f0, 'k', 'Linewidth', linewidth)
% plot(time, alt_f1, 'm', 'Linewidth', linewidth), hold off, ylim(YLim)
% xlabel('Time (sec)'), ylabel('Altitude (cm)')
% legend('baroAltitudeRaw (Debug 2)', ...
%        'baroAltitudeRaw -> 9 point median (offline)', ...
%        'baroAltitudeRaw -> pt2 filter (offline)', 'location', 'best')
% xlim(XLim)
% 
% %%
% 
% Nest     = round(5.0/Ts_baro);
% koverlap = 0.95;
% Noverlap = round(koverlap*Nest);
% window   = hann(Nest);
% delta    = 0;
% [G1, C1, freq] = my_tfestimate(alt_raw, alt_f0, window, Noverlap, Nest, Ts_baro, delta);
% [G2, C2]       = my_tfestimate(alt_raw, alt_f1, window, Noverlap, Nest, Ts_baro, delta);
% 
% G1m = tf(1,1); set(G1m, 'OutputDelay', 4.0*Ts_baro);
% G2m = tf(1,1); set(G2m, 'OutputDelay', Ts_baro);
% 
% figure(4)
% bode(G1, G2, G1m, G2m, 2*pi*freq(freq > 0 & freq <= 1/2/Ts_baro)), grid on
% legend('baroAltitudeRaw -> 9 point median (offline)', ...
%        'baroAltitudeRaw -> pt2 filter (offline)', ...
%        'Model 1 Tt = 4.0*Ts', ...
%        'Model 2 Tt = 3.0*Ts', 'location', 'best')
% 
% %%
% 
% [Alt_m, f] = my_fft(alt_f0, Ts_baro);
% [Alt_f, f] = my_fft(alt_f1, Ts_baro);
% 
% figure(5)
% plot(f, [Alt_m, Alt_f]), grid on, xlabel('Frequency (Hz)')
% set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log')
% xlim([0 1/2/Ts_baro])
