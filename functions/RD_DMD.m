function K_RD_DMD = RD_DMD(Xp,Xf,LAMBDA, N, Tfp, T1, T2)

% Computing the sparse distributed Koopman operator online as the data trickles
phi_T1     = cell(N,1);
varphi_T1  = cell(N,1);
Kp_T1      = cell(N,1);
Kp_T2      = cell(N,1);
inv_phi_T1 = cell(N,1);
inv_phi_T2 = cell(N,1);
varphi_T2  = cell(N,1);

% Initialization for Recursive part in RSD-DMD computation: 
for k1 = 1:N    
    
    % Random initialization 
    phi_T1{k1}     = 0.01*eye(size(Tfp{k1},2)); 
    varphi_T1{k1}  = zeros(size(Tfp{k1}));
    inv_phi_T1{k1} = pinv(phi_T1{k1}+LAMBDA*eye(size(phi_T1{k1})));
    
    % Koopman row blocks computation
    Kp_T1{k1}      = varphi_T1{k1}*inv_phi_T1{k1};
end

% RSD-DMD at the starting time-point
k3 = 1;
K_RD_DMD{k3} = cell2mat(Kp_T1);
% So far precomputed phi and varphi for T1 time-steps. From T1+1 onwards,
% the Koopman would be computed in real-time.

for k2 = T1:T2
    
    for k1 = 1:N % for every agent
        % Collecting new time-point
        Xp_T2 = Xp(:,k2);
        Xf_T2 = Tfp{k1}*Xf(:,k2);
        % Find phi and varphi for every new time observation
        inv_phi_T2{k1}    = inv_phi_T1{k1} - (inv_phi_T1{k1}*(Xp_T2*Xp_T2')*inv_phi_T1{k1})/(1+Xp_T2'*inv_phi_T1{k1}*Xp_T2);
        varphi_T2{k1}     = varphi_T1{k1} + Xf_T2*Xp_T2';
        Kp_T2{k1}         = varphi_T2{k1}*inv_phi_T2{k1};
        % Update phi and varphi at every agent with the current time
        % observation data which will be used in the next iteration
        inv_phi_T1{k1} = inv_phi_T2{k1};
        varphi_T1{k1}  = varphi_T2{k1};
    end
    k3 = k3 + 1;
    K_RD_DMD{k3} = cell2mat(Kp_T2);
    
end

end

