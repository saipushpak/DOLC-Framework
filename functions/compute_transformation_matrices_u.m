function [Tpp, Tfp, Tup, Rpp, Rfp, Rup] = compute_transformation_matrices_u(Adj, N, n, n_inputs)

% Inputs 
% N - Number of agents
% n - Number of states for each agent 
% Adj - Adjacency matrix of the multi-agent system 
% n_inputs - number of inputs applied at each agent 

% Outputs
% Tpp
% Tfp
% Rpp 

% Code to compute the transformation matrices required for sparse or sparse
% distributed Koopman operator

I    = eye(N); 

Tpp  = cell(N,1);
Tfp  = cell(N,1);
Tup  = cell(N,1);
Rpp  = cell(N,1);
Rfp  = cell(N,1);
Rup  = cell(N,1);

for agent_k = 1:N % for every agent
    % Define the transformation matrices for every agent
    ep = I(:,agent_k);
    
    ee = cell(N,N);
    
    for k1 = 1:N
        for k2 = 1:N
            ee{k1,k2} = zeros(n,n);
        end
    end
    
    for k1 = 1:N
        ee{k1,k1} = kron(ep(k1), eye(n));
    end
    
    Tfp{agent_k} = cell2mat(ee);
    Rfp{agent_k} = Tfp{agent_k}(any(Tfp{agent_k},2),:);
    
    %     Tfp = zeros(n,n*N);
    %     Tfp(:,2*agent_k-1:2*agent_k) = eye(n);
    
    ap = Adj(:,agent_k);
    ep = zeros(N,1);
    ep(agent_k) =1;
    
    aep = ap+ep;
    ae = cell(N,N);
    for k1 = 1:N
        for k2 = 1:N
            ae{k1,k2} = zeros(n,n);
        end
    end
    for k1 = 1:N
        ae{k1,k1} = kron(aep(k1),eye(n));
    end
    Tpp{agent_k} = cell2mat(ae);
    Rpp{agent_k} = Tpp{agent_k}(any(Tpp{agent_k},2),:);
    
    eu = cell(N,N);
    for k1 = 1:N
        for k2 = 1:N
            eu{k1,k2} = zeros(n_inputs(k1),n_inputs(k2));
        end
    end
    
    for k1 = 1:N
        eu{k1,k1} = kron(ep(k1), eye(n_inputs(k1)));
    end
    
    Tup{agent_k} = cell2mat(eu);
    Rup{agent_k} = Tup{agent_k}(any(Tup{agent_k},2),:);
end

end

