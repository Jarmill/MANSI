function out = ANSI_pair(In)
% Atomic Norm System ID
% Later to be used for Multiresolution exploration
% includes away steps to speed up convergence and promote sparsity

%% Process Input
T = In.T;
y = In.ym;
tau = In.tau.tauAtom;
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

%error_term = T*h-y;

time_elapsed = 0; 
counter = 0 ; 
tic;
if visualize || visualize_end
    figure
    dual_gap_list = [];
    dual_gap_dot_list = [];
    l1_norm = [];
end

%preprocessing
p = p_in;
A = pole_matrix(p, N, 0);

%perform scaling
scale = pole_scales(p, N)';
A_s = bsxfun(@times, A, scale');

%B = T*A; %to use the L1 formulation, everything is T*A rather than A
B = T*A_s;

%scaling needs to be done down here on the B matrix
% scale = 1./sqrt(sum(B.^2, 1))';
% B_s = bsxfun(@times, B, scale');

%pole coefficients
x = zeros([size(B, 2), 1]);
x_new = x;

pureFW = 0;

%This is just a big fancy lasso problem all things considered
%min 1/2*|T x-y|_2^2 such that |x|_A < tau_A
%turned into
%min 1/2*|T A c - y|_2^2 such that |c|_1 < tau_1
%x = A c

%here, x is used since that is pretty default


%kernel matrix 
%for fast computation of gradient
kernel = B'*B;
kernel = kernel + delta*eye(size(kernel)); %Elastic Net
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
    
    %dual_gap_dot = (-w_fw'*grad)/(norm(w_fw)*norm(grad));
    
    %% Away Step
    c = norm(x, 1);
    
    if ((t > 1) || nnz(x) > 0)
        %modify pairwise step selection if on the boundary
        
        %normal selection would go outside boundary, result in alpha_max =
        %0
        if norm(x, 1) >= tau && (sign(x(j_fw)) == at_fw || x(j_fw) == 0)
            %check these lines
            I_active_b = I_active;
            I_active_b(I_active_b == j_fw) = [];
            
            
            at_valid = -sign(x(I_active_b)); %feasible atoms, since tight on boundary
            [~, i_aw] = max(grad(I_active_b).*at_valid);
            j_aw = I_active_b(i_aw);
            at_aw = at_valid(i_aw);

        %normal pairwise step selection
        else
            [~, i_aw] = max(abs(grad(I_active)));
            j_aw = I_active(i_aw); %indices of indices. hate it.
            at_aw = -sign(grad(j_aw));
        end
    else
        j_aw = j_fw;
        at_aw = -sign(grad(j_aw));
    end
      
    w_pw = zeros([size(A, 2), 1]);
    w_pw(j_fw) = w_pw(j_fw) + at_fw*tau;
    w_pw(j_aw) = w_pw(j_aw) + at_aw*tau;

    dual_gap_fw = -w_fw'*grad;
    dual_gap_pw = -w_pw'*grad;
    
    %% Stepsize Calculation    
    %is_pair = (dual_gap_pw >= dual_gap_fw) && (norm(x, 1) < tau) && ~pureFW;
    is_pair = (dual_gap_pw >= dual_gap_fw) && ~pureFW;
    
    if is_pair     
        w = w_pw;
        if at_aw == 0
            %make sure that w has a norm of 2*tau
            w = 2*w;
        end
        
        c = norm(x, 1);
        if sign(x(j_fw)) == -sign(w(j_fw))
            c = c - 2*abs(x(j_fw));
        end

        if (sign(x(j_aw)) == -sign(w(j_aw)))
            c = c - 2*abs(x(j_aw));
        end

        alpha_max = (tau - c)/(2*tau);

    else
        w = w_fw;
        alpha_max = 1;
    end

    dual_gap = -w'*grad;
    dual_gap_dot = dual_gap / (norm(w, 2)*norm(grad, 2));

    
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
    
    
    %terminate = dual_gap_dot < 0.01;
    %terminate = dual_gap_dot < 0.001;
    %terminate = dual_gap_dot < 1e-5;
    %terminate = dual_gap < 1e-5;
    terminate = dual_gap < 1e-8 || (t == t_max-1);
    
    if visualize || visualize_end
        %record data for later
        dual_gap_list(t) = dual_gap;
        dual_gap_dot_list(t) = dual_gap_dot;
        l1_norm(t) = norm(x_new, 1);
    end
    
    if visualize || (visualize_end && terminate)
        %plots and stuff        
        clf
        subplot(4, 3, [1, 2, 4, 5])
        hold on
        active_old = find(x ~= 0);
        active_new = find(x_new ~= 0);
        stem3(real(p(active_old)), imag(p(active_old)), x(active_old), '.')
        stem3(real(p(active_new)), imag(p(active_new)), x_new(active_new), '.')
        
        %active stuff
        scatter(real(p), imag(p), [], abs(grad), '.')

        cb = colorbar;
        ylabel(cb, 'abs(grad)')
        
        best_pole = p(j_fw);
        if is_pair
            best_pole_pair = p(j_aw);
        else
            best_pole_pair = NaN;
        end

        scatter3(real(best_pole), imag(best_pole), x_new(j_fw), [], 'r')     
        if is_pair
            scatter3(real(best_pole_pair), imag(best_pole_pair), x_new(j_aw), [], 'xg')
        end
        
        th = linspace(0, 2*pi, 400);
        plot3(cos(th), sin(th), zeros(size(th)), 'color', [0 .5 0]);
        
        view(3)
        title(strcat('Atom Weights(', num2str(nnz(x_new)), '/', num2str(length(x_new)), ')'))
        xlabel('Re(z)')
        ylabel('Im(z)')
        zlabel('Residue')
        hold off
        
        
        
        %optimal stepsize (optional?)
        subplot(4, 3, [3, 6])
        %box on
        alpha_N = 301;
        %alpha_list = linspace(0, 1, alpha_N);
        alpha_list = logspace(-5, 0, alpha_N);
        alpha_pts = x + w*alpha_list;
        alpha_val = sum((B*alpha_pts - y).^2, 1) /2;
        
        alpha_best_val = norm(B*(x + alpha*w) - y)^2 /2;
        
        hold on
        plot(alpha_list, alpha_val);
        scatter(alpha, alpha_best_val, '*')
        hold off
        title(strcat('alpha stepsize (\alpha=', num2str(alpha, 3), ')'))
        xlabel('alpha')
        ylabel('f(x+alpha*w)')
        set(gca, 'xscale', 'log')
        
        
        %duality gap plot
        subplot(4, 3, 9)
        semilogy(dual_gap_list)
        title('Duality Gap (f(x)-f(x*) <= DG')                
        %semilogy(dual_gap_dot_list)       
        %title('Duality Gap (Forward) Cos(Angle)')
        xlabel('iterations')
        ylabel('Duality Gap')
        
        %l1_norm
        subplot(4, 3, 12)
        hold on
        plot(l1_norm)
        plot([1, t], [tau, tau], '--k')
        hold off
        ylabel('|Sx|')
        title('L1 norm')
        hold off
        
        %impulse response
        subplot(4, 3, [7, 8, 10, 11])
        hold on
        y_hat = B*x_new;
        
    plot(y,'*');
    hold on;
    plot(y_hat,'o');
    hold off;
    legend groundtruth estimated
    xlabel('t');
    ylabel('system response');
    title('Atomic norm approximation');
    
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

x_out = scale .* x;

%final output
out.c = x_out;
%impulse response
%could probably speed this up
h_out = A*x_out;
y_out = T*h_out;
out.h = h_out;
out.y = y_out;

end
