function [C] = MAVSC(n, K, mu,lambda1)

num_iter = 50;
max_mu = 1e6;
rho = 1.5;
err_thr = 10^-4;
num_views = 3;

C1 = repmat({zeros(n,n)}, 1, num_views);
C2 = repmat({zeros(n,n)}, 1, num_views);
A = repmat({zeros(n,n)}, 1, num_views);
A_prev = repmat({zeros(n,n)}, 1, num_views);
Y1 = repmat({zeros(n,n)}, 1, num_views);
Y2 = repmat({zeros(n,n)}, 1, num_views);
C_centroid = zeros(n,n);

alpha = rand(1,num_views);
alpha = alpha/sum(alpha,2);

mu1 = mu;
mu2 = mu;
mu = [mu1 mu2];

iter = 0;
converged = false;
sumalpha=0;

while iter < num_iter && ~converged
    
    for v = 1:num_views
        A_prev{v} = A{v}; % save previous value
        [C1{v}, C2{v}, Y1{v}, Y2{v}, A{v}] = oneview...
            (n, K{v}, C1{v}, C2{v}, C_centroid, Y1{v}, Y2{v}, ...
            alpha(1,v), lambda1, mu);
    end
    % update num_views
    for v = 1:num_views
        alpha(1,v) = 0.5/norm(C_centroid-C1{v},'fro');
    end
    
    % update centroid
    for v = 1:num_views
        C_centroid = C_centroid+alpha(1,v)*C1{v};
    end
    for v = 1:num_views
        sumalpha = sumalpha + alpha(1,v);
    end
    C_centroid = C_centroid/sumalpha;

    % check convergence
    converged = true;
    for v = 1:num_views
        err1 = max(max(abs(A_prev{v}-A{v})));
        err2 = max(max(abs(A{v}-C1{v}+diag(diag(C1{v})))));
        err3 = max(max(abs(A{v}-C2{v})));

       if err1>err_thr || err2>err_thr || err3>err_thr
           converged = false;
           break
       end    
    end
    iter = iter + 1;
    mu = min(rho*mu,max_mu);
end

C = C_centroid;
end
