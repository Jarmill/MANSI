%coefficients
c = [1 2 -1 3];
Fs = 1;

%poles
pa = [0.8+0.4j, 1.0j, -0.5, 0.8-0.4j]; 

pc = conj(pa);

[pu, ip, ipc] = union(pa, pc);

order = length(pu);
cu = zeros([1, length(c)]);

%find where values in p are in the union
%this is hackery
[~, ipsort] = sort(pa);

p_locations = find(ismember(pu, pa));
cu(p_locations) = c(ipsort);

%fill in values of the residues
r = zeros(size(pu));
r(imag(pu) == 0) = cu(imag(pu) == 0);

ind_trig = find(imag(pu) ~= 0);
ind_cos = ind_trig(1:2:length(ind_trig));
ind_sin = ind_trig(2:2:length(ind_trig));

residue_top    = (cu(ind_cos) + 1.0j*cu(ind_sin))/2;
residue_bottom = conj(residue_top);

r(ind_cos) = residue_top;
r(ind_sin) = residue_bottom;

%transfer function testing
% [b, a] = residue(r, pu, 0);
% b = real(b);
% sys = tf(b, a, 1);
[z, k] = zeros_from_poles(r, pu);

sys = zpk(z, pu, k, Fs);

figure
pzmap(sys)