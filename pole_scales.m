function [ scales ] = pole_scales(p, N, complex)
%scaling coefficients for the set of poles p
%this is a nontrivial task. used to normalize pole response

%N is the horizon. If N>=0, this is the number 
%If N < 0, use infinite horizon.

if nargin < 3
    complex = 0;
end

scales = zeros(size(p));

r_all = abs(p);
r2_all = r_all.^2;
theta_all = angle(p);

%'horizon length'. is a very tricky quantity.
base = N-1;
%base = N+1;
if complex
    norm2_complex = zeros(size(p));
    index   = find(r_all < 1);
    index_b = find(r_all >= 1);
    
    complex_numerator = 1;
    complex_denominator = 1 - r2_all(index);
    
    if N > 0 && N ~= Inf
        %finite horizon
        complex_numerator = complex_numerator - r2_all(index).^(base);
    end
    
    %pole is inside unit circle
    norm2_complex(index) = complex_numerator./complex_denominator;
    
    %pole is on boundary of unit circle
    norm2_complex(index_b) = base;
    
    scales = 1./sqrt(norm2_complex);
        
else

    %find the values of the sums;
    %all poles
    index_exp = find(imag(p) == 0 & r_all < 1);
    index_cos  = find(imag(p) > 0 & r_all < 1);
    index_sin  = find(imag(p) < 0 & r_all < 1);

    %index_inside = [index_exp; index_cos; index_sin];
    index_inside = [index_exp, index_cos, index_sin];
    num_inside = length(index_inside);
    r2 = r2_all(index_inside);
    theta = theta_all(index_inside);

    %special care needs to be taken about poles on unit circle, L'hopital's
    %rule.

    index_exp_b = find(imag(p) == 0 & r_all == 1);
    index_cos_b = find(imag(p) >  0 & r_all == 1);
    index_sin_b = find(imag(p) <  0 & r_all == 1);

    index_boundary = [index_cos_b, index_sin_b];
    theta_b = theta_all(index_boundary);

    %all things together
    index_exp_all = [index_exp, index_exp_b];
    index_cos_all = [index_cos, index_cos_b];
    index_sin_all = [index_sin, index_sin_b];


    %% Finding the values of the sums r^(2n) and r^(2n)*cos(2n theta)
    %power sum r^(2n)
    ps_numerator   = ones(size(index_inside));
    ps_denominator = 1 - r2;
    
    %trig sum r^(2n)*cos(2n theta)
    ts_numerator   = 1 - r2 .* cos(2*theta);
    ts_denominator = 1 - 2* r2 .* cos(2*theta) + r2.^2;
    
    %trig sum on boundary (limit as r->1)
    ts_numerator_b = 1-cos(2*theta_b);
    ts_denominator_b = 2*(1 - cos(2*theta_b));

    if N > 0 && N ~= Inf
        %finite horizon
%         ps_numerator = ps_numerator - r2.^(N+1);
        %ts_numerator = ts_numerator + r2.^(N+1).*(-cos(2*(N+1)*theta) + r2 .* cos(2*N*theta));

        ps_numerator = ps_numerator - r2.^(base);
        ts_numerator = ts_numerator + r2.^(base).*(-cos(2*(base)*theta) + r2 .* cos(2*(base-1)*theta));

%         ts_numerator = ts_numerator + r2.^(base).*(-cos(2*(base)*theta) + r2 .* cos(2*(base-1)*theta));
        
        %when r->1 and this are sinusoids
        %needs some tlc
        %ts_numerator_b = ts_numerator_b + cos(2*(N)*theta_b) - cos(2*(N+1)*theta_b);
        ts_numerator_b = ts_numerator_b + cos(2*(base-1)*theta_b) - cos(2*(base)*theta_b);

    end

    
    power_sum_i = ps_numerator./ps_denominator;
    trig_sum_i = ts_numerator./ts_denominator;

    %power_sum_b = N+1;
    power_sum_b = base;
    trig_sum_b = ts_numerator_b ./ ts_denominator_b;

    %thrown values into the arrays
    power_sum = zeros(size(p));
    trig_sum = zeros(size(p));

    power_sum(index_inside) = power_sum_i;
    trig_sum(index_inside)  = trig_sum_i;

    power_sum(index_boundary) = power_sum_b;
    trig_sum(index_boundary)  = trig_sum_b;

    power_sum(index_exp_b) = base;
    trig_sum(index_exp_b)  = 0;


    %% Using sums to evaluate scaling factors

    %L2 norm squared of pole responses
    norm2_exp = power_sum(index_exp_all);
    norm2_cos = (power_sum(index_cos_all) + trig_sum(index_cos_all))./2;
    norm2_sin = (power_sum(index_sin_all) - trig_sum(index_sin_all))./2;

    % Sum of r^(2n)
    scales(index_exp_all) = 1./sqrt(norm2_exp);

    % Sum of (r^(n)*cos(n theta))^2
    scales(index_cos_all) = 1./sqrt(norm2_cos);

    % Sum of (r^(n)*sin(n theta))^2
    scales(index_sin_all) = 1./sqrt(norm2_sin);
end

end

