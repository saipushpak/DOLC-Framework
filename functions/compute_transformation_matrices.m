function [Tpp, Tfp, Rpp] = compute_transformation_matrices(Adj, N, n)

% Inputs 
% N - Number of agents
% n - Number of states for each agent 
% Adj - Adjacency matrix of the multi-agent system 

% Outputs
% Tpp
% Tfp
% Rpp 

% Code to compute the transformation matrices required for sparse or sparse
% distributed Koopman operator

Tpp  = cell(N,1);
Tfp  = cell(N,1);
Rpp  = cell(N,1);

for agent_k = 1:N % for every agent
    
    % Define the transformation matrices for every agent
    Tfp{agent_k} = zeros(n,n*N);
    Tfp{agent_k}(:,2*agent_k-1:2*agent_k) = eye(n);
    
    ap = Adj(:,agent_k);
    ep = zeros(N,1);
    ep(agent_k) =1;
    
    aep = ap+ep;
    ae = cell(N,N);
    for k2 = 1:N
        for k3 = 1:N
            ae{k2,k3} = zeros(n);
        end
    end
    for k2 = 1:N
        ae{k2,k2} = kron(aep(k2),eye(n));
    end
    Tpp{agent_k} = cell2mat(ae);
    Rpp{agent_k} = Tpp{agent_k}(any(Tpp{agent_k},2),:);
end


end

