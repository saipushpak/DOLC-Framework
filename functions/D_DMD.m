function K_D_DMD = D_DMD(Xp,Xf,LAMBDA, N, Tfp)
% D-DMD: Robust Distributed DMD 

% Inputs 
% Xp and Xf are one time-step separated datasets
% Xp - 1 --> t-1
% Xf - 2 --> t
% LAMBDA - regularization factor 
% N - Number of agents (assuming each node as an agent) 
% Tfp - Transformation matrix to collect Xf dataset corresponding to each
% agent 

% Output: Koopman operator using Robust distributed DMD

% Identify the Koopman in a distributed manner 
% Can be solved analytically

% Collect time-series data corresponding to each agent
X_fp = cell(N,1);
K_p  = cell(N,1);
jk   = zeros(N,1);

for k1 = 1:N % for every agent
    X_fp{k1,1} = Tfp{k1}*Xf;       
    
    tic
    Yf = X_fp{k1,1}*Xp';
    Yp = Xp*Xp';
    
    K_p{k1,1} = (Yf)*pinv(Yp + LAMBDA*eye(size(Yp)));
    jk(k1) = toc;
end

K_D_DMD = cell2mat(K_p);

end

