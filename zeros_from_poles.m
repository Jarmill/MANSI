function [z, k] = zeros_from_poles(r, p)
%find the zeros of the transfer function given the set of residues and
%poles
%H(z) = Prod(z-si)/Prod(z-pi) = Sum ri/(z-pi). Find si, roots of polynomial
%b.
%all poles are unique, no double poles. may be close together though,
%nothing to be done about that.
%a more stable/applied version of residue.m
N = length(p);

b = zeros([1, length(p)]);

if N > 1
    for idx = 1:N
        p_active = p(idx);
        ptemp = p;
        ptemp(idx) = NaN; %omit the pole at idx, partial fraction
        btemp = poly(ptemp);
        b = b + r(idx) .* btemp;
    end
else
    b = r;
end

%b is now the polynomial on top
%find the roots of b to get the zeroes

%b might still have complex parts, on order of 1e-16. kill them.
b = real(b);

%might not have a full set of zeroes, which is ok.
first_nonzero = find(b~=0, 1, 'first');
b = b(first_nonzero:end);

%find the gain of the system, and rescale
k = b(1);
b = b / k;

%break down system into roots
z = roots(b);
%z(imag(z) < 1e-15) = real(z(imag(z) < 1e-15)); %more numerical errors

end