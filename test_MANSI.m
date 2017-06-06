% system identification of trajectories

addpath(genpath('../3rdParty/thesisCode'));

% rng(0);
% data_generation;
% y = data(2001, :)';
% roots(sys_par{1})

% load dict;
% p = rp.*exp(1i*thetap);
% In.p_in = p;

% p = genPoleOnRing(0.1, 0.95, 0.1, 0.1);
% In.p_in = p;
% In.k = length(p);

% p = [0.3+0.5j, 0.3-0.5j, 0.8];
% p = [0.5*exp(0.1j), 0.5*exp(-0.1j), 0.8];
p = [0.5];
%p = [-0.005 + 0.5j; -0.005 - 0.5j];
%p = [0.95j; -0.95j];
%p = [-0.5 + 0.5j, -0.5 - 0.5j, 0.7];
b = [1];
a = poly(p);
sysd = tf(b, a, 1);
[c_true, ~, ~] = residue(b, a);
y = impulse(sysd, 0:100);

N = size(y, 1);
In.ym = y;
In.T = eye(N);
%In.tau.tauAtom = 5;
In.tau.tauAtom = 2;
%In.tau.tauAtom = 1;
In.t_max = 1000;
In.k = 150;


%radius = 8;
%radius = 11;
Npoles = 2*radius - 1;
[poles_xx, poles_yy] = meshgrid(linspace(-1, 1, Npoles));
poles = poles_xx + 1.0j*poles_yy;
poles_circ = poles(abs(poles) <= 1);
In.p_in = poles_circ';
In.k = length(poles_circ);

In.visualize = 1;
%Out = atomic_SISO(In);
%Out = ANSI_forward(In);
Out = ANSI_away(In);
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
stem3(real(p), imag(p), c_true)
stem3(real(poles_active), imag(poles_active), c(active_ind))

th = linspace(0, 2*pi, 400);
plot3(cos(th), sin(th), zeros(size(th)), 'color', [0 .5 0] )

hold off
axis square
view(3)
title('Pole Map')
%legend groundtruth estimated
xlabel('Re(z)')
ylabel('Im(z)')
zlabel('Coefficients')
% [~, indind] = sort(-abs(c));
% c = c(indind);
% p = Out.p;
% p = p(indind);
% 
% ind = find(c);
% p_est = p(ind);