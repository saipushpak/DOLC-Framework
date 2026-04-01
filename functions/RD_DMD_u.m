function [A, B, learning_error] = RD_DMD_u(Xp, Xf, U, LAMBDA, N, Tfp, Rfp, T1, T2, n, which_initialization)

% This code assumes there is control that can be applied at every agent

% Computing the distributed Koopman operator online as the data trickles
phi_T1          = cell(N,1);
varphi_T1       = cell(N,1);
Ap_T1           = cell(N,1);
Ap_T2           = cell(N,1);
inv_phi_T1      = cell(N,1);
inv_phi_T2      = cell(N,1);
varphi_T2       = cell(N,1);
Bp_T1           = cell(N,1);
Bp_T2           = cell(N,1);

% Initialization for Recursive part in RSD-DMD computation:
for k1 = 1:N
    if strcmp(which_initialization, 'random')
        % Random initialization
        phi_T1{k1}     = 0.01*eye(n*N+N); % [x u]*[x u]' number of states per agent + total number of states
    elseif strcmp(which_initialization, 'initial_window')
        Z_k = [Xp(:,1:T1);U(:,1:T1)];
        phi_T1{k1}     = Z_k*Z_k';
    end
    varphi_T1{k1}  = zeros(n+n*N,n*N+N); % [y x]*[x u]'
    inv_phi_T1{k1} = pinv(phi_T1{k1}+LAMBDA*eye(size(phi_T1{k1})));

    % Koopman row blocks computation
    temp              = varphi_T1{k1}*inv_phi_T1{k1};
    Ap_T1{k1}      = temp(1:n,1:n*N);
    Bp_T1{k1}      = temp(1:n,1+n*N:end);

end

% RSD-DMD at the starting time-point
A     = cell(T2-T1,1);
B     = cell(T2-T1,1);
k3    = 1;
A{k3} = cell2mat(Ap_T1);
B{k3} = cell2mat(Bp_T1);
% So far precomputed phi and varphi for T1 time-steps. From T1+1 onwards,
% the Koopman would be computed in real-time.

learning_error.A = NaN(T2-T1,1);
learning_error.B = NaN(T2-T1,1);

for k2 = T1:T2

    for k1 = 1:N % for every agent
        % Collecting new time-point
        Xp_T2 = [Xp(:,k2); U(:,k2)];
        Xf_T2 = [Rfp{k1}*Tfp{k1}*Xf(:,k2); Xp(:,k2)];
        % Find phi and varphi for every new time observation
        inv_phi_T2{k1}    = inv_phi_T1{k1} - (inv_phi_T1{k1}*(Xp_T2*Xp_T2')*inv_phi_T1{k1})/(1+Xp_T2'*inv_phi_T1{k1}*Xp_T2);
        varphi_T2{k1}     = varphi_T1{k1} + Xf_T2*Xp_T2';
        temp                  = varphi_T2{k1}*inv_phi_T2{k1};
        Ap_T2{k1}          = temp(1:n,1:n*N);
        Bp_T2{k1}          = temp(1:n,1+n*N:end);
        % Update phi and varphi at every agent with the current time
        % observation data which will be used in the next iteration
        inv_phi_T1{k1} = inv_phi_T2{k1};
        varphi_T1{k1}  = varphi_T2{k1};
    end
    k3    = k3 + 1;
    A{k3} = cell2mat(Ap_T2);
    B{k3} = cell2mat(Bp_T2);

    learning_error.A(k2) = norm(A{k3-1} - A{k3}, 'fro');
    learning_error.B(k2) = norm(B{k3-1} - B{k3}, 'fro');

end

end

