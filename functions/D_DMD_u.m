function [A, B] = D_DMD_u(Xp, Xf, U, LAMBDA, N, Tfp, Rfp, n)

% This code assumes there is control that can be applied at every agent 

Ap = cell(N,1);
Bp = cell(N,1);

for agent_k = 1:N % for every agent
    
    X_bar = [Xp; U];
    Y_bar = [Rfp{agent_k}*Tfp{agent_k}*Xf; Xp];
    
    V = Y_bar*X_bar';
    G = X_bar*X_bar';
    
    M = V*pinv(G + LAMBDA*blkdiag(eye(n*N), zeros(N)));
    
    Ap{agent_k} = M(1:n, 1:n*N);
    Bp{agent_k} = M(1:n, n*N+1:end);
        
end

A = cell2mat(Ap);
B = cell2mat(Bp);

end

