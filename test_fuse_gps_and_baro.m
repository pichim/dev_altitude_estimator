clc, clear all
%%

is_matlab = false;  % false -> octave

if is_matlab
  addpath ../fcn_bib
else
  addpath fcn_bib_copy
  pkg load control
end

Ts = 1/40;
downsample_divider = 4;

% create a similar function like ctzsnooze
N = 25*downsample_divider;
time = (0:N-1).' *Ts;
alt = -cos(2*pi*0.3*time) + 0.3;
alt = cumtrapz(time, alt);
alt = alt - alt(1);
alt = [zeros(10,1); alt(1:N-10)]*10;
dalt_end = mean(diff(alt(N-10:N))) / Ts;
alt = [alt; alt(end) + Ts*dalt_end + time(1:N/4) * dalt_end];
N = N + N/4;
time = (0:N-1).' *Ts;
Gf = c2d(tf(1, [1/(2*pi*1) 1]), Ts, 'tustin');
alt = filter(Gf.num{1}, Gf.den{1}, alt);

% create sensor readings baro_alt and gps_alt
baro_alt = alt;
gps_alt  = zeros(N,1);
downsample_cntr = 0;
for i = 1:N
    if mod(downsample_cntr, downsample_divider) == 0
        altitude_gps_i = alt(i,:);
    end
    gps_alt(i,:) = altitude_gps_i;
    downsample_cntr = downsample_cntr + 1;
end

% add  noise and/or offset
apply_corruption = 1;
baro_offset = apply_corruption * ones(N,1);
baro_noise = apply_corruption * cumtrapz(time, sqrt(0.1) * randn(N,1));
baro_alt = baro_alt + baro_offset + baro_noise;
gps_noise = apply_corruption * cumtrapz(time, sqrt(0.1) * randn(N,1));
gps_alt = gps_alt + gps_noise;

% fuse the data accoring to the trust values, this will be the dc
gps_trust = 0.7;
alt_dc_fused = gps_trust * gps_alt + (1 - gps_trust) * baro_alt;

% fuse the dc with a complementary filter and baro
f_cut = 0.3; % Hz
Gpt3 = tf(get_filter('pt3', f_cut, Ts));
alt_fused = filter(Gpt3.num{1}, Gpt3.den{1}, alt_dc_fused) + baro_alt - filter(Gpt3.num{1}, Gpt3.den{1}, baro_alt);

lw = 1.2;

figure(1)
plot(time, alt, 'LineWidth', lw), grid on, hold on
stairs(time, [baro_alt, gps_alt, alt_dc_fused], 'LineWidth', lw), hold off
xlabel('Time (sec)'), ylabel('Altitude (m)'), xlim([0 time(end)])
legend('ideal altitude', 'baro alt', 'gps alt', 'alt dc fused', 'location', 'northwest')

figure(2)
plot(time, alt, 'LineWidth', lw), grid on, hold on
stairs(time, [alt_dc_fused, alt_fused], 'LineWidth', lw), hold off
xlabel('Time (sec)'), ylabel('Altitude (m)'), xlim([0 time(end)])
legend('ideal altitude', 'alt dc fused', 'alt fused', 'location', 'northwest')

figure(3)
stairs(time(1:end-1), diff([alt, alt_dc_fused, alt_fused]), 'LineWidth', lw), grid on
xlabel('Time (sec)'), ylabel('d/dt Altitude (m)'), xlim([0 time(end)])
legend('ideal altitude', 'alt dc fused', 'alt fused', 'location', 'northwest')
