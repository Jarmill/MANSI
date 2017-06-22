function [ A ] = pole_matrix( p, N, complex)
%Vandermonde matrix on the poles
%performs task of finding impulse responses components to horizon

%p: set of poles in complex plane (must include conjugates)
%N: horizon, number of samples
if nargin < 3
    complex = 0;
end

k = length(p);

%response of 1/(z-p) = r^n (cos(n theta) + j sin(n theta))
r = abs(p)';
theta = angle(p)';

%scale = (1 - r.^2);
% r2 = r.^2;
% recip_scale = (1 - r2.^(N)) ./ (1 - r2);
% scale = sqrt(1./(recip_scale));

%raise r to powers
index_N = ones(1,N-2);
%mag = [ones(size(r)); cumprod(r(index_N, :))];
mag = (r.^(0:N-2))';
%angles = (0:N-2)'*theta; %n*theta
angles = (0:N-2)'*theta'; %n*theta

%optional step
angles = mod(angles, 2*pi);

%poles with a+bi will have cos, a-bi will have sin
%exponentials are normal, with cos(0) = 1
%angles(:, theta < 0) = angles(:, theta < 0) + pi/2;

if complex
    trig_angles = exp(1.0j*angles);
else
    trig_angles = zeros(size(angles));
    trig_angles(:, imag(p) >= 0) = cos(angles(:, imag(p) >= 0));
    trig_angles(:, imag(p) < 0) = sin(angles(:, imag(p) < 0));
end
trig_angles(abs(trig_angles) < 1e-15) = 0;

vandermonde = mag .* trig_angles;

%A = [zeros(1,k);ones(1,k); vandermonde] * diag(scale);
A = [zeros(1,k); vandermonde];


end

