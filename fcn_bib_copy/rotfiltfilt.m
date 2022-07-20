% rotating filtfilt-function
function xf = rotfiltfilt(Bf, Af, sinarg, x)

% signal size
[Nx, nx] = size(x);
xf = zeros(Nx, nx);
p = exp(1j*(sinarg));

for i = 1:nx
    % eliminate mean
    y = x(:,i) - mean(x(:,i));
    yR = y.*p;
    yQ = y.*conj(p);
    % filtering in transformed coordinates
    yR = filtfilt(Bf, Af, yR);
    yQ = filtfilt(Bf, Af, yQ);
    % back transformation
    % xf(:,i) = real((yR.*conj(p) + yQ.*p)*0.5);
    xf(:,i) = real((yR.*conj(p) + yQ.*p));
end

end