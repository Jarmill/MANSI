% system identification of trajectories

addpath(genpath('../3rdParty/thesisCode'));

% rng(0);
% data_generation;
% y = data(2001, :)';
% roots(sys_par{1})

% load dict;
% p_sys = rp.*exp(1i*thetap);
% In.p_in = p_sys;

% p_sys = genPoleOnRing(0.1, 0.95, 0.1, 0.1);
% In.p_in = p_sys;
% In.k = length(p_sys);

% p_sys = [0.3+0.5j, 0.3-0.5j, 0.8];
% p_sys = [0.5*exp(0.1j), 0.5*exp(-0.1j), 0.8];
%p_sys = [0.5*exp(0.1j), 0.5*exp(-0.1j)];
%p_sys = [-0.9; 0.5];
%p_sys = [0; 0.9];
%p_sys = [0.65];
%p_sys = [1/sqrt(2)];
%p_sys = [-0.005 + 0.5j; -0.005 - 0.5j];
%p_sys = [0.95j; -0.95j];
%p_sys = [0.7j; -0.7j; 0.2];
%p_sys = [-0.5 + 0.5j, -0.5 - 0.5j];
%p_sys = [-1];
%p_sys = [0.75; 0.95];
%p_sys = [0.0; 0.8];
%p_sys = [0.3; 0.5];
%p_sys = [0.3; 0.7];
p_sys = [0.6];
%p_sys = [-0.5 + 0.5j, -0.5 - 0.5j, 0.7];

b = [1];
%b = [length(p_sys) sum(p_sys)]
%b = [1 0.5];
%b = [1 0];
a = poly(p_sys);
Fs = 1;
sysd = tf(b, a, Fs);
[r, ~, ~] = residue(b, a);
c_true = r;
%deal with complex poles/complex residues
%as per convention, imag(p)>0: cos, imag(p)<0: sin, imag(p)=0: exp
% c_true(imag(p_sys) > 0) = 2*real(c_true(imag(p_sys) > 0));
% c_true(imag(p_sys) < 0) = 2*imag(c_true(imag(p_sys) < 0));

c_true(imag(c_true) > 0) = 2*real(c_true(imag(c_true) > 0));
c_true(imag(c_true) < 0) = 2*imag(c_true(imag(c_true) < 0));

In.visualize = 1;
In.visualize_end = 1;

%In.tau.tauAtom = 2.25;
%In.tau.tauAtom = 5;
%In.tau.tauAtom = 1.2;
In.tau.tauAtom = 1.6;
%In.tau.tauAtom = 1;


In.tau.delta = 1e-4; %Elastic Net Regularization
%In.tau.delta = 0; %Elastic Net Regularization

In.tau.lambda = 1e-1;
%In.tau.lambda = 1e-2;

In.t_max = 500;
In.k = 150;

N = 101;
response_type = 0;
if response_type == 1
    %step response
    %no idea how this will go
    
    %the scaling for this is screwed up. massively.
    %and the closed-form formulas are horrific
    In.T = tril(ones(N));
    y = step(sysd, 0:N-1);
else
    %impulse response
    In.T = eye(N);
    y = impulse(sysd, 0:N-1);
end
In.ym = y;

%radius = 5;
radius = 10;
%radius = 30;
%radius = 40;
%radius = 200;
Npoles = 2*radius + 1;
rho = 1;

%poles on unit circle will ring forever, are excluded.
[poles_xx, poles_yy] = meshgrid(linspace(-rho, rho, Npoles));
poles = poles_xx + 1.0j*poles_yy;
poles_circ = poles(abs(poles) <= 1);
%poles_circ = reshape(poles_circ, [1, length(poles_circ)]);

%real axis only, exponents for testing
%poles_circ = linspace(-1, 1, Npoles)';

%clear up numerical artifacts
poles_circ(abs(imag(poles_circ)) < 1e-15) = real(poles_circ(abs(imag(poles_circ)) < 1e-15));
poles_circ(abs(real(poles_circ)) < 1e-15) = 1.0j * imag(poles_circ(abs(real(poles_circ)) < 1e-15));

In.p_in = poles_circ';
In.k = length(poles_circ);


%Out = atomic_SISO(In);
%Out = ANSI_forward(In);
%Out = ANSI_away(In);
%Out = ANSI_pair(In);
Out = ADMMSI(In);
y_hat = Out.h;

% In.Nf = 2^8;
% In.idx_avail = 1:N;
% Out = atomic_SISO_FFT_xikang(In);
% y_hat = Out.h;

figure
subplot(1, 2, 1)
plot(y,'*');
hold on;
plot(y_hat,'o');
hold off;
legend groundtruth estimated
xlabel('t');
ylabel('input and output signal');
title('Atomic norm approximation');

c = Out.c;
subplot(1, 2, 2)
active_ind = find(c ~= 0);
poles_active = poles_circ(active_ind);

hold on
stem3(real(p_sys), imag(p_sys), c_true)
stem3(real(poles_active), imag(poles_active), c(active_ind))

th = linspace(0, 2*pi, 400);
plot3(cos(th), sin(th), zeros(size(th)), 'color', [0 .5 0] )

hold off
axis square
view(3)
title(strcat('Pole Map (', num2str(nnz(c)), '/', num2str(length(c)),')'))

%legend groundtruth estimated
xlabel('Re(z)')
ylabel('Im(z)')
zlabel('Coefficients')
