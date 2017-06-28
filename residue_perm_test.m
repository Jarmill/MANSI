%coefficients
Fs = 1;

%poles
pa = [0.8+0.4j, 1.0j, -0.5, 0.8-0.4j]; 
c = [1 2 -1 3];

% pa = [0.8 + 0.4j, 1.0j, -0.5];
% c = [1 -3 2];

% pa = [-0.8+.2j, -0.8-.2j];
% c =  [2, -1];

%pa = [-0.8j, 0.2, -0.4, 0.8j];
%c =  [2, -3, 0.5, 1];
 
% pa = [-0.8j, -0.2j, 0.8j];
% c =  [1, 2, -3];

pc = conj(pa);

[pu, ip, ipc] = union(pa, pc);

order = length(pu);
cu = zeros([1, length(pu)]);

%find where values in p are in the union
%this is hackery
[~, ipsort] = sort(pa);

p_locations = find(ismember(pu, pa));
cu(p_locations) = c(ipsort);

%fill in values of the residues
r = zeros(size(pu));
r(imag(pu) == 0) = cu(imag(pu) == 0);

ind_real = find(imag(pu) == 0);
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

%try combining complex conjugates
sys_accum = 0;

num_real = length(ind_real);
num_trig = length(ind_trig)/2;

%real poles
for i = ind_real
    sys_curr = zpk([], pu(i), r(i), Fs);
    sys_accum = sys_accum + sys_curr;
end

%complex poles
for i = ind_cos
    rc = r(i);
    pcurr = pu(i);
    k = rc + conj(rc);
    s = (rc*conj(pcurr) + conj(rc)*pcurr);
    
    %k = 0, pure cos
    %k = Inf, pure sin
    %k = else, mix
    %need to validate
    
    if k == 0
        sys_curr = zpk([], [pcurr, conj(pcurr)], -s, Fs);
    else
        sys_curr = zpk(s/k, [pcurr, conj(pcurr)], k, Fs);
    end

    sys_accum = sys_accum + sys_curr;
end

figure
subplot(2, 1, 1)
pzmap(sys)
title('original sys')

subplot(2, 1, 2)
pzmap(sys_accum)
title('accumulated sys')