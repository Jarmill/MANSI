function [ A ] = pole_matrix_complex( p, N)
%Vandermonde matrix on the poles
%performs task of finding impulse responses components to horizon
%complex version, will have complex numbers

%p: set of poles in complex plane (must include conjugates)
%N: horizon, number of samples

%raise r to powers
% index_N = ones(1,N-2);
% mag = [ones(size(r)); cumprod(r(index_N, :))];
% 
% vandermonde = mag .* trig_angles;
% 
% A = [zeros(1,k); vandermonde];

k = length(p);


index_N = ones(1, N-2);

mag = cumprod(p(index_N, :));

%A = mag;
A = [zeros(1,k); ones(1,k); mag];

end

