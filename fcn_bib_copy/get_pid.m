function [Cpi, Cd] = get_pid(PID, Ts)

Kp = PID(1); Ki = PID(2); Kd = PID(3); 
Cpi = ss(Kp + Ki*Ts*tf([1 0], [1 -1], Ts));
Cd  = ss(Kd/Ts*tf([1 -1], [1 0], Ts));

end

