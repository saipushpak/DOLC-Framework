function [A, B, C] = DMD_u(Xp,Xf,U,LAMBDA,N,n)
% Robust DMD: Robust Dynamic Mode Decomposition with control 

% Inputs 
% Xp and Xf are one  time-step separated datasets
% Xp - 1 --> t-1
% Xf - 2 --> t
% U  - 1 --> t-1 (control input datad) 
% LAMBDA - regularization factor 
% N - Number of agents (assuming each node as an agent) 
% n - number of states at each agent 

% Output:
% A - Koopman operator (system matrix)
% B - input matrix
% C - output matrix (for DMD case, it must be identity matrix)

% Form V and G matrices and compute A, B, C matrices
X_bar = [Xp; U]; 
Y_bar = [Xf; Xp]; 

V = Y_bar*X_bar';
G = X_bar*X_bar'; 

M = V*pinv(G + LAMBDA*eye(size(G)));
A = M(1:n*N,1:n*N);
B = M(1:n*N,n*N+1:end);
C = M(n*N+1:end,1:n*N);

end