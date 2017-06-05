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
sysd = tf([1 ], poly(p), 1);
y = impulse(sysd, 0:100);

N = size(y, 1);
In.ym = y;
In.T = eye(N);
%In.tau.tauAtom = 5;
In.tau.tauAtom = 1;
In.t_max = 1000;
In.k = 150;


radius = 8;
Npoles = 2*radius - 1;
[poles_xx, poles_yy] = meshgrid(linspace(-1, 1, Npoles));
poles = poles_xx + 1.0j*poles_yy;
poles_circ = poles(abs(poles) <= 1);
In.p_in = poles_circ';
In.k = length(poles_circ);

In.visualize = 1;
%Out = atomic_SISO(In);
Out = ANSI_forward(In);
y_hat = Out.h;

% In.Nf = 2^8;
% In.idx_avail = 1:N;
% Out = atomic_SISO_FFT_xikang(In);
% y_hat = Out.h;


plot(y,'*');
hold on;
plot(y_hat,'o');
hold off;
legend groundtruth estimated
xlabel('t');
ylabel('input and output signal');
title('Atomic norm approximation');


% c = Out.c;
% [~, indind] = sort(-abs(c));
% c = c(indind);
% p = Out.p;
% p = p(indind);
% 
% ind = find(c);
% p_est = p(ind);