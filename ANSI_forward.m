function out = ANSI_forward(In)
% Atomic Norm System ID
% Later to be used for Multiresolution exploration

%% Process Input
T = In.T;
y = In.ym;
tau = In.tau.tauAtom;
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

if isfield(In, 'visualize') && In.visualize == 1
    visualize = 1;
else
    visualize = 0;
end

%error_term = T*h-y;

time_elapsed = 0; 
counter = 0 ; 
tic;
if visualize
    figure
    dual_gap_dot_list = [];
    l1_norm = [];
end

%preprocessing
p = p_in;
A = pole_matrix(p, N);
B = T*A; %to use the L1 formulation, everything is T*A rather than A

%pole coefficients
x = zeros([size(B, 2), 1]);
x_new = x;

%This is just a big fancy lasso problem all things considered
%min 1/2*|T x-y|_2^2 such that |x|_A < tau_A
%turned into
%min 1/2*|T A c - y|_2^2 such that |c|_1 < tau_1
%x = A c

%here, x is used since that is pretty default


%kernel matrix 
%for fast computation of gradient
kernel = B'*B;
Bty = B'*y;

%main loop of processing
for t = 1:t_max
    %% calculate gradient
    %bottleneck of computation
    I_active = find(abs(x) > 0);
    x_a = x(I_active); %nonzero x values

    %the current bottleneck of the code
    %still can't figure out a faster formulation
    %grad = kernel(:, I_active)*x_a - Atb;
    grad = kernel(:, I_active)*x_a;
    grad = grad - Bty;
    
    %% Forward Step
    [~, j_fw] = max(abs(grad));
    if(length(j_fw) > 1)
        j_fw = j_fw(randi(1, length(j_fw)));
    end
    
    w_fw = -x;
    at_fw = -sign(grad(j_fw));
    w_fw(j_fw) = w_fw(j_fw) + at_fw*tau;
    
    dual_gap = -w_fw'*grad;
    dual_gap_dot = (-w_fw'*grad)/(norm(w_fw)*norm(grad));
    
    %% Stepsize Calculation
    %forward only
        is_away = 0;
        w = w_fw;
        alpha_max = 1;
        j_curr = j_fw;
        
    
    I_active_w = find(abs(w) > 0); %union is inefficient apparently
    w_a = w(I_active_w); %nonzero w values
    
    alpha_top = w_a'*grad(I_active_w);
    alpha_bottom = sum(w_a .* (kernel(I_active_w, I_active_w)*w_a)); %w'A'Aw
        if alpha_bottom == 0
        break
    end
    
    alpha = -alpha_top/alpha_bottom;
    
    %make sure alpha is between 0 and 1
    alpha = max(0, min(alpha_max, alpha)); 
    
    %% Reporting
    x_new = x + alpha * w;
    
    terminate = dual_gap_dot < 0.01;
    
    if visualize
        %things happen here
        dual_gap_dot_list(t) = dual_gap_dot;
        l1_norm(t) = norm(x_new, 1);
        
        clf
        subplot(5, 1, [1, 2, 3])
        hold on
        stem3(real(p), imag(p), x)
        stem3(real(p), imag(p), x_new)
        
        scatter(real(p), imag(p), [], grad, 'x')
        colorbar
        
        p_best = p(j_curr);
        scatter(real(p_best), imag(p_best), [], 'r')     
        
        th = linspace(0, 2*pi, 400);
        plot3(cos(th), sin(th), zeros(size(th)));
        
        view(3)
        hold off
        
        subplot(5, 1, 4)
            semilogy(dual_gap_dot_list)
        title('Duality Gap (Forward) Cos(Angle)')
        xlabel('iterations')
        ylabel('Angle between -w and grad')
        
        subplot(5, 1, 5)
        hold on
        plot(l1_norm)
        plot([1, t], [tau, tau], '--k')
        hold off
        ylabel('|x|')
        title('L1 norm')
    end
    
    %% Update
    x = x_new;
    if terminate
        break
    end
    
end

%timing
time_elapsed = toc;
out.time_elapsed = time_elapsed;
out.iter = t;

%final output
out.c = x;
%impulse response
%could probably speed this up
h_out = A*x;
y_out = T*h_out;
out.h = h_out;
out.y = y_out;

end
