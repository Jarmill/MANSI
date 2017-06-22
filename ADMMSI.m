function out = ADMMSI(In)
% ADMM System ID
% I know that ADMM works on the elastic net, so this will be a test.
% pole responses over a grid can be preprocessed, which should massively
% speed up the full-rank update in the least squares prox operator.

%% Process Input
T = In.T;
b = In.ym;
%tau = In.tau.tauAtom;
lambda = In.tau.lambda;
delta = In.tau.delta;
t_max = In.t_max;

m = In.k;

N = size(T,2);

if isfield(In,'h0')
    h = In.h0;
else
    h = zeros(N,1);
end

if isfield(In,'p_in')
    p_in = In.p_in;
    pole_array = true;
else
    pole_array =false;
end

if isfield(In,'complex')
    if In.complex
        complex = true;
    else
        complex = false;
    end
else
    complex =false;
end

if isfield(In, 'visualize') && In.visualize == 1
    visualize = 1;
else
    visualize = 0;
end

if isfield(In, 'visualize_end') && In.visualize_end == 1
    visualize_end = 1;
else
    visualize_end = 0;
end


time_elapsed = 0; 
tic;
%visualization prep
if visualize || visualize_end
    figure
    comp_slack_list = [];
    en_error = [];
    gap_error = [];
    
    th = linspace(0, 2*pi, 400);
    th_cos = cos(th);
    th_sin = sin(th);
    th_z   = zeros(size(th));
end

%% preprocessing
p = p_in;
A = pole_matrix(p, N, complex);

%perform scaling
scale = pole_scales(p, N, complex);
%scale_true = 1./sqrt(sum(A.*conj(A), 1)); %for reference
A_s = bsxfun(@times, A, scale);

%B = T*A; %to use the L1 formulation, everything is T*A rather than A
B = T*A_s;

%kernel matrix 
kernel = B'*B;
kernel = kernel + delta*eye(size(kernel)); %Elastic Net
Btb = B'*b;

%admm setup
t = 0;
t_max = 200;

rho = 1e-5;
rho_max = 10000;
rho_increment = 1.5;

%variables
x = zeros([size(A, 2), 1]);
y = zeros(size(x));
z = zeros(size(x));

%errors and termination
old_en_error = Inf;
old_gap_error = Inf;
terminate = 0;

%% main loop of processing
while ~terminate
    %least squares update on elastic net kernel
    %can be sped up by preprocessed eigendecomposition
    x = prox_least_sq(z - y/rho, kernel, Btb, rho);
    
    %l1 update
    z = prox_l1(x + y/rho, lambda/rho);
    
    %lagrange update
    y = y + rho*(x - z);
    
    %error value tracking
    new_en_error = elastic_net_error(B, b, lambda, delta, (x+z)/2);
    %new_en_error = elastic_net_lagrangian_error(B, b, lambda, delta, x, z, y);
    rel_en_error = abs(new_en_error - old_en_error)/old_en_error;
    
    new_gap_error = norm(x-z, 2)^2;
    rel_gap_error = abs(new_gap_error - old_gap_error)/old_gap_error;
    
    %check for ending
    %terminate = (t > t_max) || ((rel_en_error < 3e-5) && (new_gap_error < 1e-10));
    terminate = (t > t_max) || (new_gap_error < 1e-11);
    %visualization
    if visualize || (terminate && visualize_end)
        clf
        %visualization goes here
        %poles (x/z)
        subplot(3, 3, [1, 2, 4, 5])
        hold on
        p_active_x = find(x ~= 0);
        p_active_z = find(z ~= 0);
        stem3(real(p(p_active_x)), imag(p(p_active_x)), x(p_active_x), '.')
        stem3(real(p(p_active_z)), imag(p(p_active_z)), z(p_active_z), '.')
    
        %lagrange multipliers
        scatter(real(p), imag(p), [], y, '.')
        cb = colorbar;
        ylabel(cb, 'y')
        
        %unit circle
        plot3(th_cos, th_sin, th_z, 'color', [0 .5 0]);
        
        view(3)
        title(strcat('pole atom weights (\rho=', num2str(rho, '%10.2e'), ', z:', num2str(nnz(z)), '/', num2str(length(z)), ')'))
        xlabel('Re(z)')
        ylabel('Im(z)')
        zlabel('Coefficient')
        hold off
        
        
        %complementary slackness
        comp_slack = abs(y'*(x-z));
        comp_slack_list(t+1) = comp_slack;
        subplot(3, 3, 3)
        plot(comp_slack_list)
        xlabel('iterations')
        ylabel('<y,(x-z)>')
        title('|Complementary Slackness|')
        set(gca, 'yscale', 'log')
        
        %errors
        error(t+1) = new_en_error;
        gap_error(t+1) = new_gap_error;
        
        subplot(3, 3, 6)
        plot(error)
        xlabel('iterations')
        ylabel('elastic net error')
        set(gca,'yscale','log')
        
        yyaxis right
        plot(gap_error)
        ylabel('consensus disagreement')
        set(gca,'yscale','log')
        title('elastic net error')
        
        %responses
        subplot(3, 3, [7, 8, 9])
        hold on
        plot(b, '*')
        plot(B*x, 'o')
        plot(B*z, 'x')
        xlabel('t')
        ylabel('system response')
        title('ADMM system response')
        legend('true', 'x (non-sparse)', 'z (sparse)', 'Location', 'northeast')
        hold off
    end

    %setup for next iteration
    t = t + 1;
    rho = min(rho*rho_increment, rho_max);
    old_en_error = new_en_error;
    old_gap_error = new_gap_error;
end

%% Output
%timing
time_elapsed = toc;
out.time_elapsed = time_elapsed;
out.iter = t;

%x is non-sparse, while z is always sparse. use z as the output.

x_out = scale' .* z;

%final output
out.c = x_out;
%impulse response
%could probably speed this up
h_out = A*x_out;
y_out = T*h_out;
out.h = h_out;
out.y = y_out;

end
