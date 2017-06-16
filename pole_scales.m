function [ scales ] = pole_scales( p, N )
%scaling coefficients for the set of poles p
%this is a nontrivial task. used to normalize pole response

%N is the horizon. If N>=0, this is the number 
%If N < 0, use infinite horizon.
scales = zeros(size(p));

%find the values of the sums;
%all poles
r_all = abs(p);
r2_all = r_all.^2;
theta_all = angle(p);


index_exp = find(imag(p) == 0 & r_all < 1);
index_cos  = find(imag(p) > 0 & r_all < 1);
index_sin  = find(imag(p) < 0 & r_all < 1);

index_inside = [index_exp index_cos index_sin];
num_inside = length(index_inside);
r2 = r2_all(index_inside);
theta = theta_all(index_inside);

%special care needs to be taken about poles on unit circle, L'hopital's
%rule.

index_exp_b = find(imag(p) == 0 & r_all == 1);
index_cos_b = find(imag(p) >  0 & r_all == 1);
index_sin_b = find(imag(p) <  0 & r_all == 1);

index_boundary = [index_cos_b index_sin_b];
theta_b = theta_all(index_boundary);

%all things together
index_exp_all = [index_exp index_exp_b];
index_cos_all = [index_cos index_cos_b];
index_sin_all = [index_sin index_sin_b];


%% Finding the values of the sums r^(2n) and r^(2n)*cos(2n theta)

ps_numerator   = ones(1, length(index_inside));
ts_numerator   = 1 - r2 .* cos(2*theta);

ts_numerator_b = 1-cos(2*theta_b);

if N > 0
    %finite horizon
    ps_numerator = ps_numerator - r2.^(N+1);
    ts_numerator = ts_numerator + r2.^(N+1).*(-cos(2*(N+1)*theta) + r2 .* cos(2*N*theta));
    
    %when r->1 and this are sinusoids
    ts_numerator_b = ts_numerator_b + cos(2*N*theta_b) - cos(2*(N+1)*theta_b);
    
end

ps_denominator = 1 - r2;
ts_denominator = 1 - 2* r2 .* cos(2*theta) + r2.^2;

ts_denominator_b = 2*(1 - cos(2*theta_b));

power_sum_i = ps_numerator./ps_denominator;
trig_sum_i = ts_numerator./ts_denominator;

power_sum_b = N+1;
trig_sum_b = ts_numerator_b ./ ts_denominator_b;

%thrown values into the arrays
power_sum = zeros(size(p));
trig_sum = zeros(size(p));

power_sum(index_inside) = power_sum_i;
trig_sum(index_inside)  = trig_sum_i;

power_sum(index_boundary) = power_sum_b;
trig_sum(index_boundary)  = trig_sum_b;

power_sum(index_exp_b) = N+1;
trig_sum(index_exp_b)  = 0;


%% Using sums to evaluate scaling factors

%L2 norm squared of pole responses
norm2_exp = power_sum(index_exp_all);
norm2_cos = (power_sum(index_cos_all) + trig_sum(index_cos_all))/2;
norm2_sin = (power_sum(index_sin_all) - trig_sum(index_sin_all))/2;

% Sum of r^(2n)
scales(index_exp_all) = 1./sqrt(norm2_exp);

% Sum of (r^(n)*cos(n theta))^2
scales(index_cos_all) = 1./sqrt(norm2_cos);

% Sum of (r^(n)*sin(n theta))^2
scales(index_sin_all) = 1./sqrt(norm2_sin);


end

