function [Cpi, Cd, Gf, PID, para_used] = get_tfs_from_para(para, ind_ax, Ts_cntr, throttle_avg)

if nargin == 3
    throttle_avg = 0;
end
dterm_lpf1_dyn_expo = 10;
gyro_lpf1_dyn_expo = 10;

% param mapping, example: 20210827_flight_00.txt
% logging rate, sampling time, 
%       "looptime",125                  gyro:    Ts     = 125 mus            <-> 8 kHz
%       "pid_process_denom",2           cntrl:   Ts_cnt = 2*Ts     = 250 mus <-> 4 kHz
%       "frameIntervalPDenom",2         logging: Ts_log = 2*Ts_cnt = 500 mus <-> 2 kHz 
% const. gyro lowpass 1
%       "gyro_soft_type",0              type: 0:PT1 1:BIQUAD    
%       "gyro_lowpass_hz",310           cutoff frequency: 0:OFF
% const. gyro lowpass 2
%       "gyro_soft2_type",0             type: 0:PT1 1:BIQUAD
%       "gyro_lowpass2_hz",310          cutoff frequency: 0:OFF
% const. gyro notch fitler 1 and 2
%       "gyro_notch_hz","111,135"       centerFreq: 0:off
%       "gyro_notch_cutoff","102,120"   cutoffFreq: 0:off
% const. dterm lowpass 1
%       "dterm_filter_type",0           type: 0:PT1 1:BIQUAD    
%       "dterm_lpf_hz",209              cutoff frequency: 0:OFF
% const. dterm lowpass 2
%       "dterm_filter2_type",0          type: 0:PT1 1:BIQUAD
%       "dterm_lpf2_hz",209             cutoff frequency: 0:OFF
% const. dterm notch fitler
%       "dterm_notch_hz",93             centerFreq: 0:off
%       "dterm_notch_cutoff",86         cutoffFreq: 0:off
% iterm relax
%       "iterm_relax",1                 0:off 1:on
% controller parameters
%       "rollPID","53,64,43,0"
%       "pitchPID","65,54,57,0"
%       "yawPID","58,58,0,0"

% Gf: y -> yf: gyro filter
Gf = ss(tf(1, 1, Ts_cntr));
filter_enumeration = {'pt1', 'biquad', 'pt2', 'pt3'};
% gyro lowpass filter 1
if para.gyro_lowpass_hz > 0
    para_used.gyro_lowpass_hz = para.gyro_lowpass_hz;
    para_used.gyro_soft_type  = para.gyro_soft_type;
    Gf = Gf * get_filter(filter_enumeration{para.gyro_soft_type + 1}, ...
                         para.gyro_lowpass_hz, ...
                         Ts_cntr);
end
% dynamic gyro lowpass filter 1
if para.gyro_lowpass_dyn_hz(1) > 0
    % make sure Gf is 1 at start, this is not possible in current bf
    Gf = ss(tf(1, 1, Ts_cntr));
    para_used.gyro_lowpass_dyn_hz = para.gyro_lowpass_dyn_hz;
    para_used.gyro_soft_type      = para.gyro_soft_type;
    para_used.gyro_lpf_hz_avg     = get_fcut_from_exp(para.gyro_lowpass_dyn_hz(1), ...
                                                      para.gyro_lowpass_dyn_hz(2), ...
                                                      gyro_lpf1_dyn_expo, ...
                                                      throttle_avg);
    para_used.dterm_lpf_throttle_avg = throttle_avg;
    Gf = Gf * get_filter(filter_enumeration{para.dterm_filter_type + 1}, ...
                         para_used.dterm_lpf_hz_avg, ...
                         Ts_cntr);
end
% gyro lowpass filter 2
if para.gyro_lowpass2_hz > 0
    para_used.gyro_lowpass2_hz = para.gyro_lowpass2_hz;
    para_used.gyro_soft2_type  = para.gyro_soft2_type;
    Gf = Gf * get_filter(filter_enumeration{para.gyro_soft2_type + 1}, ...
                         para.gyro_lowpass2_hz, ...
                         Ts_cntr);
end
% gyro notch filter 1
if para.gyro_notch_hz(1) > 0
    para_used.gyro_notch_hz(1)     = para.gyro_notch_hz(1);
    para_used.gyro_notch_cutoff(1) = para.gyro_notch_cutoff(1);
    Gf = Gf * get_filter('notch', ...
                         [para.gyro_notch_cutoff(1), para.gyro_notch_hz(1)], ...
                         Ts_cntr);
end
% gyro notch filter 2
if para.gyro_notch_hz(2) > 0
    para_used.gyro_notch_hz(2)     = para.gyro_notch_hz(2);
    para_used.gyro_notch_cutoff(2) = para.gyro_notch_cutoff(2);
    Gf = Gf * get_filter('notch', ...
                         [para.gyro_notch_cutoff(2), para.gyro_notch_hz(2)], ...
                         Ts_cntr);
end

% Gd: d/dt(yf) -> d/dt(yf)f: dterm filter
Gd = ss(tf(1, 1, Ts_cntr));
% filter_enumeration = {'pt1', 'biquad', 'pt2', 'pt3'};
% dterm lowpass filter 1
if para.dterm_lpf_hz > 0
    para_used.dterm_lpf_hz      = para.dterm_lpf_hz;
    para_used.dterm_filter_type = para.dterm_filter_type;
    Gd = Gd * get_filter(filter_enumeration{para.dterm_filter_type + 1}, ...
                         para.dterm_lpf_hz, ...
                         Ts_cntr);
end
% dynamic dterm lowpass filter 1
if para.dterm_lpf_dyn_hz(1) > 0
    % make sure Gd is 1 at start, this is not possible in current bf
    Gd = ss(tf(1, 1, Ts_cntr));
    para_used.dterm_lpf_dyn_hz  = para.dterm_lpf_dyn_hz;
    para_used.dterm_filter_type = para.dterm_filter_type;
    para_used.dterm_lpf_hz_avg  = get_fcut_from_exp(para.dterm_lpf_dyn_hz(1), ...
                                                    para.dterm_lpf_dyn_hz(2), ...
                                                    dterm_lpf1_dyn_expo, ...
                                                    throttle_avg);
    para_used.dterm_lpf_throttle_avg = throttle_avg;
    Gd = Gd * get_filter(filter_enumeration{para.dterm_filter_type + 1}, ...
                         para_used.dterm_lpf_hz_avg, ...
                         Ts_cntr);
end
% dterm lowpass filter 2
if para.dterm_lpf2_hz > 0
    para_used.dterm_lpf2_hz      = para.dterm_lpf2_hz;
    para_used.dterm_filter2_type = para.dterm_filter2_type;
    Gd = Gd * get_filter(filter_enumeration{para.dterm_filter2_type + 1}, ...
                         para.dterm_lpf2_hz, ...
                         Ts_cntr);
end
% dterm notch filter
if para.dterm_notch_hz > 0
    para_used.dterm_notch_hz     = para.dterm_notch_hz;
    para_used.dterm_notch_cutoff = para.dterm_notch_cutoff;
    Gd = Gd * get_filter('notch', ...
                         [para.dterm_notch_cutoff, para.dterm_notch_hz], ...
                         Ts_cntr);
end

% PID parameters, Todo: make this not so messy
switch ind_ax
    case 1
        if size(para.rollPID, 2) == 5
            para.rollPID = para.rollPID([1 2 3 5]);
        end 
        if size(para.rollPID, 2) == 4
            PID = para.rollPID  .* [0.032029, 0.244381    , 0.000529, 0];
        else
            PID = para.rollPID  .* [0.032029, 0.244381    , 0.000529];
        end
    case 2
        if size(para.pitchPID, 2) == 5
            para.pitchPID = para.pitchPID([1 2 3 5]);
        end 
        if size(para.pitchPID, 2) == 4
            PID = para.pitchPID .* [0.032029, 0.244381    , 0.000529, 0];
        else
            PID = para.pitchPID .* [0.032029, 0.244381    , 0.000529];
        end
    case 3
        if size(para.yawPID, 2) == 5
            para.yawPID = para.yawPID([1 2 3 5]);
        end 
        if size(para.yawPID, 2) == 4
            PID = para.yawPID   .* [0.032029, 0.244381*2.5, 0.000529, 0];
        else
            PID = para.yawPID   .* [0.032029, 0.244381*2.5, 0.000529];
        end
    otherwise
end
[Cpi, Cd] = get_pid(PID, Ts_cntr);
Cd = Cd * Gd;

end

