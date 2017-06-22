N = 101;
t = 0:N-1;
radius = 4;
%radius = 20;
%radius = 30;
%radius = 40;
%radius = 200;
Npoles = 2*radius + 1;

%rho = 0.95;
rho = 1;

%poles on unit circle will ring forever, are excluded.
[poles_xx, poles_yy] = meshgrid(linspace(-rho, rho, Npoles));
poles = poles_xx + 1.0j*poles_yy;
poles_circ = poles(abs(poles) <= 1);

%real axis only, exponents for testing
%poles_circ = linspace(-1, 1, Npoles)';

%clear up numerical artifacts
poles_circ(abs(imag(poles_circ)) < 1e-15) = real(poles_circ(abs(imag(poles_circ)) < 1e-15));
poles_circ(abs(real(poles_circ)) < 1e-15) = 1.0j * imag(poles_circ(abs(real(poles_circ)) < 1e-15));

A = pole_matrix_complex(poles_circ', N);
scale = pole_scales(poles_circ, N, 1);

A_s = bsxfun(@times, A, scale');

figure
plot3(t, real(A_s), imag(A_s))
xlabel('t')
ylabel('Re(\rho^n)')
zlabel('Im(\rho^n)')