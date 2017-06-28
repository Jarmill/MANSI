function [sys, sys_reduce] = zpk_from_poles(c, p, tol, Fs)
%from coefficients and poles, generates system
%needs to find residues of expansion first 

%c: list of coefficients of responses (real, sin/cos)
%p: locations of active poles in system
%tol: tolerance of pole/zero cancellation
%Fs: sample rate of system

ind_p_real = find(abs(imag(p)) <= 1e-8);
ind_p_trig = find(abs(imag(p))  > 1e-8);

N_real = length(ind_p_real);
p_real = p(ind_p_real);
r_real = c(ind_p_real);

p_trig = p(ind_p_trig);
p_trigc = conj(p_trig);

%full set of complex poles
ptu = union(p_trig, p_trigc);

%here's where things get odd
%find poles with nonzero coefficients
%load these in to the c array 
%generate residues r from c
c_trig = zeros(size(ptu));
N_trig = length(ptu);
[~, ipsort] = sort(p_trig);
p_locations = find(ismember(ptu, p_trig));
c_trig(p_locations) = c(ipsort);

%fill in values of residues
r_trig = zeros(size(ptu));

ind_cos = 1:2:N_trig;
ind_sin = 2:2:N_trig;

residue_top    = (c_trig(ind_cos) + 1.0j*c_trig(ind_sin))/2;
residue_bottom = (c_trig(ind_cos) - 1.0j*c_trig(ind_sin))/2;

r_trig(ind_cos) = residue_top;
r_trig(ind_sin) = residue_bottom;

%% find resultant system
sys = 0;
sys_reduce = 0;

for i = 1:N_real
    rc = r_real(i);
    pc = p_real(i);
    sys_curr_real = zpk([], pc, rc, Fs);
    sys = sys + sys_curr_real;
    sys_reduce = minreal(sys_reduce + sys_curr_real, tol);
end

for i = 1:2:N_trig
    rc = r_trig(i);
    pc = ptu(i);
    k = rc + conj(rc);
    s = (rc*conj(pc) + conj(rc)*pc);
    
    %k = 0, pure cos
    %k = Inf, pure sin
    %k = else, mix
    %need to validate
    
    if k == 0
        sys_curr_trig = zpk([], [pc, conj(pc)], -s, Fs);
    else
        sys_curr_trig = zpk(s/k, [pc, conj(pc)], k, Fs);
    end
    
    sys = sys + sys_curr_trig;
    sys_reduce = minreal(sys_reduce + sys_curr_trig, tol);
end

end