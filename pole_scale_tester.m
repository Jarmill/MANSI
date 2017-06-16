%p  =  [0.5000 + 0.5000i   0.5000 - 0.5000i   0.3000 + 0.0000i];
p = [0.0322 + 0.9990i]; %really close to jw axis/unit circle

%N = 5;
N = 101;
scale = 0;
A = pole_matrix(p, N, scale);
l2A = sqrt(sum(A.^2, 1));
spA = 1./l2A;

sp  = pole_scales(p, N-1);
%sp  = pole_scales(p, -1);
%sp  = pole_scales(p, N);
% spN = pole_scales(p, N);

new_norm_sp = sp .* l2A
% new_norm_spN = spN .* l2A