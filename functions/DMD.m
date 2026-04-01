function K_DMD = DMD(Xp,Xf,LAMBDA)
% Robust DMD: Robust Dynamic Mode Decomposition

% Inputs 
% Xp and Xf are one  time-step separated datasets
% Xp - 1 --> t-1
% Xf - 2 --> t
% LAMBDA - regularization factor 

% Output: Koopman operator using Robust DMD
Yf = Xf*Xp';
Yp = Xp*Xp';

K_DMD = (Yf)*pinv(Yp + LAMBDA*eye(size(Yp)));
end