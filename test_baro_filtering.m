clc, clear all
%%

try
  addpath ../fcn_bib
catch
  addpath fcn_bib_copy
end

load data_20220807 % save data_20220807 data time

% ind_pid_p   =  3: 5;
% ind_pid_i   =  6: 8;
% ind_pid_d   =  9:10;
% ind_pid_f   = 11:13;
% ind_rc      = 14:16;
ind_setp    = 18:21;
% ind_baro    = 24;
% ind_gyro_f  = 26:28;
% ind_acc     = 29:31;
% ind_mot     = 36:39;
ind_debug   = 32:35;
% ind_pid_sum = 48:50;
% ind_pid_err = 55:57;
% ind_pi_sum  = 62:64;

% DEBUG_SET(DEBUG_BARO, 0, state);
% DEBUG_SET(DEBUG_BARO, 1, filteredBaroAltitudeRaw - baroGroundAltitude);
% DEBUG_SET(DEBUG_BARO, 2, baroAltitudeRaw - baroGroundAltitude);
% DEBUG_SET(DEBUG_BARO, 3, baro.BaroAlt); // metres * 100
% 0: 
% 1: raw baro to alt conversion with median and pt3
% 2: raw baro to alt conversion
% 3: old baro value

figure(1)
subplot(211)
stairs(time(1:end-1), diff(time))
subplot(212)
stairs(time, data(:, ind_debug(1)))

mean(diff(time))
median(diff(time))

%%

Ts_log = 1/2e3;

% Ndata = size(data, 1);
% ind_eval = false(Ndata,1);
% time_0 = 0;
% for i = 1:Ndata
%     if time(i) - time_0 > 0.02 && data(i,ind_debug(1)) > 2.5
%         ind_eval(i) = true;
%         time_0 = time(i);
%     end
% end
% figure(99)
% subplot(211)
% stairs(diff(time(ind_eval)))

Ts_baro = 1/40;
ind_eval = 20:Ts_baro/Ts_log:size(data, 1);
alt_raw  = data(ind_eval, ind_debug(3));
alt      = data(ind_eval, ind_debug(2));
alt_old  = data(ind_eval, ind_debug(4));
throttle = data(ind_eval, ind_setp(4));
time = time(ind_eval);

data_ = [time, alt_raw, alt, alt_old, throttle];
     
table_data = array2table(data_);
table_data.Properties.VariableNames = {'time (sec)', ...
    'DEBUG_BARO_2', 'DEBUG_BARO_1', 'DEBUG_BARO_3', 'Throttle', };

writetable(table_data, 'LOG00001_20220728_ctzsnooze_PR_#11775.csv')

%%

N = size(alt_raw, 1);

alt_f0 = alt_raw;
n = 9;
for i = n:N
    alt_window = alt_raw(i-n+1:i);
    [~, idx] = sort(alt_window);
    if mod(n, 2) == 0
        alt_f0(i) = mean(alt_window(idx(n/2:n/2+1)));
    else
        alt_f0(i) = alt_window(idx((n+1)/2));
    end
end
% alt_f0 = filter([0.5 0.5], [1 0], alt_f0);

alt_f1 = alt_raw;
% n = 7;
% for i = n:N
%     alt_window = alt_raw(i-n+1:i);
%     [~, idx] = sort(alt_window);
%     if mod(n, 2) == 0
%         alt_f1(i) = mean(alt_window(idx(n/2:n/2+1)));
%     else
%         alt_f1(i) = alt_window(idx((n+1)/2));
%     end
% end

f_cut = 2.0;
[G, B, A] = get_filter('pt2', f_cut, Ts_baro);
alt_f1 = filter(B, A, alt_f1);

G_fast = get_filter('pt2', f_cut, 1/120);
figure(99)
bode(G, G_fast), grid on

linewidth = 1.0;

% XLim = [63 68]; YLim = [9200 10600];
% XLim = [153 161]; YLim = [13200 14600];
% XLim = [395 400]; YLim = [1400 2600];
XLim = [25 38]; YLim = [3400 6600];
% XLim = [0 inf]; YLim = [-inf inf];

figure(2)
ax(1) = subplot(221);
stairs(time, [alt_raw, alt, alt_old], 'Linewidth', linewidth), grid on, ylim(YLim)
ylabel('Altitude (cm)')
legend('baroAltitudeRaw (Debug 2)', 'filteredBaroAltitudeRaw (Debug 1)', 'baro.BaroAlt (Debug 3)', 'location', 'best')
ax(2) = subplot(223);
stairs(time, throttle, 'k', 'Linewidth', linewidth), grid on
xlabel('Time (sec)'), ylabel('Throttle')
ax(3) = subplot(122);
stairs(time, alt_raw, 'Linewidth', linewidth), grid on, hold on
stairs(time, alt_f0, 'k', 'Linewidth', linewidth)
stairs(time, alt_f1, 'm', 'Linewidth', linewidth), hold off, ylim(YLim)
xlabel('Time (sec)'), ylabel('Altitude (cm)')
legend('baroAltitudeRaw (Debug 2)', ...
       'baroAltitudeRaw -> 9 point median (offline)', ...
       'baroAltitudeRaw -> pt2 filter (offline)', 'location', 'best')
linkaxes(ax, 'x')
xlim(XLim)

figure(3)
plot(time, alt_raw, 'Linewidth', linewidth), grid on, hold on
plot(time, alt_f0, 'k', 'Linewidth', linewidth)
plot(time, alt_f1, 'm', 'Linewidth', linewidth), hold off, ylim(YLim)
xlabel('Time (sec)'), ylabel('Altitude (cm)')
legend('baroAltitudeRaw (Debug 2)', ...
       'baroAltitudeRaw -> 9 point median (offline)', ...
       'baroAltitudeRaw -> pt2 filter (offline)', 'location', 'best')
xlim(XLim)

%%

Nest     = round(5.0/Ts_baro);
koverlap = 0.95;
Noverlap = round(koverlap*Nest);
window   = hann(Nest);
delta    = 0;
[G1, C1, freq] = my_tfestimate(alt_raw, alt_f0, window, Noverlap, Nest, Ts_baro, delta);
[G2, C2]       = my_tfestimate(alt_raw, alt_f1, window, Noverlap, Nest, Ts_baro, delta);

G1m = tf(1,1); set(G1m, 'OutputDelay', 4.0*Ts_baro);
G2m = tf(1,1); set(G2m, 'OutputDelay', Ts_baro);

figure(4)
bode(G1, G2, G1m, G2m, 2*pi*freq(freq > 0 & freq <= 1/2/Ts_baro)), grid on
legend('baroAltitudeRaw -> 9 point median (offline)', ...
       'baroAltitudeRaw -> pt2 filter (offline)', ...
       'Model 1 Tt = 4.0*Ts', ...
       'Model 2 Tt = 3.0*Ts', 'location', 'best')

%%

[Alt_m, f] = my_fft(alt_f0, Ts_baro);
[Alt_f, f] = my_fft(alt_f1, Ts_baro);

figure(5)
plot(f, [Alt_m, Alt_f]), grid on, xlabel('Frequency (Hz)')
set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log')
xlim([0 1/2/Ts_baro])
