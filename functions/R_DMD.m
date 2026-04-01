function K_R_DMD = R_DMD(Xp,Xf,LAMBDA, N, T1, T2)
% D-DMD: Recursive DMD 

% Inputs 
% Xp and Xf are one time-step separated datasets
% Xp - 1 --> t-1
% Xf - 2 --> t
% LAMBDA - regularization factor 
% N - Number of agents (assuming each node as an agent) 
% T1 - Given number of observations (doesn't really play a role, it is only
% here for the sake of completeness) 
% T2 - total number of time steps, the recursive Koopman is computed 

% Output
% K_R_DMD - cell containing T2 number of Koopman operators wrt each time point 

phi_T1      = 0.01*eye(2*N);
varphi_T1   = zeros(2*N);
inv_phi_T1  = pinv(phi_T1+LAMBDA*eye(2*N));
k3          = 1;
K_R_DMD{k3} = varphi_T1*inv_phi_T1;

for k2 = T1:T2
    k3          = k3 + 1;
    % Collecting new time-point
    Xp_T2       = Xp(:,k2);
    Xf_T2       = Xf(:,k2);
    % Find phi and varphi for every new time observation
    inv_phi_T2  = inv_phi_T1 - (inv_phi_T1*(Xp_T2*Xp_T2')*inv_phi_T1)/(1+Xp_T2'*inv_phi_T1*Xp_T2);
    varphi_T2   = varphi_T1 + Xf_T2*Xp_T2';
    K_R_DMD{k3} = varphi_T2*inv_phi_T2;
    % Update phi and varphi at every agent with the current time
    % observation data which will be used in the next iteration
    inv_phi_T1  = inv_phi_T2;
    varphi_T1   = varphi_T2;    
end

end