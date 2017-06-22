function ts = trig_sum(p, N)
%evaluates the sum r^(2n) * Cos(2n theta)
r = abs(p);
theta = angle(p);
r2 = r.^2;
%ts_numerator = 1;
%ts_denominator = 1 - r2;

ts_numerator   = 1 - r2 .* cos(2*theta);
ts_denominator = 1 - 2* r2 .* cos(2*theta) + r2.^2;

if N > 0 && N ~= Inf
    base = N-1;
%    ts_numerator = ts_numerator + r2.^(base).*(-cos(2*(base)*theta) + r2 .* cos(2*(base-1)*theta));
    ts_numerator = ts_numerator + r2.^(N-1).*(-cos(2*(N-1)*theta) + r2 .* cos(2*(N-2)*theta));

end

ts = ts_numerator ./ ts_denominator;


end