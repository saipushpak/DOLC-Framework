function K_S_D_DMD = SD_DMD(Xp,Xf,LAMBDA, N, Tpp, Tfp, Rpp)
% S-D-DMD: Robust Sparse Distributed DMD 

% Inputs 
% Xp and Xf are one time-step separated datasets
% Xp - 1 --> t-1
% Xf - 2 --> t
% LAMBDA - regularization factor 
% N - Number of agents
% Tpp 
% Tfp
% Rpp 

% Output: Koopman operator using Robust sparse distributed DMD

% Identify the Koopman in a distributive way
% Can be solved analytically
% Leverage the Adjacency matrix

% Collect time-series data corresponding to each agent
X_fp = cell(N,1);
X_pp = cell(N,1);
K_p  = cell(N,1);

jk = zeros(N,1);
for k1 = 1:N % for every agent
    X_fp{k1,1} = Tfp{k1}*Xf;
    X_pp{k1,1} = Rpp{k1}*Tpp{k1}*Xp;    
    tic
    Yf = X_fp{k1,1}*X_pp{k1,1}';
    Yp = X_pp{k1,1}*X_pp{k1,1}';
    
    K_p{k1,1} = (Yf)*pinv(Yp + LAMBDA*eye(size(Yp)))*Rpp{k1};
    jk(k1) = toc;
end

K_S_D_DMD = cell2mat(K_p);

end