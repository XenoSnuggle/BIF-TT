% This is an example usage of inchowrm BIF-TT
clear all;
param.dt = 0.2; % time step size
param.N = 20; % number of time steps
param.beta = 1.0; 
param.wc = 2.5;
param.wmax = 4* param.wc;
param.L = 400;
param.epsilon = 1.0;
param.Delta = 1.0;
param.xi = 0.8; % Kondo parameter
param.M = 1; % inchworm truncation parameter
param.svd_tol = 1e-8;
param.ram_range = [200,200]; % limit the storage size of BIF-TT
% param.ttrounding_tol = 1e-6;

method = 'tol'; % or 'rank' : the choice of the way of TT-rounding

Hs = zeros(2,2);
Hs(1,1) = param.epsilon; Hs(2,2) = -param.epsilon;
Hs(1,2) = param.Delta; Hs(2,1) = param.Delta;
bare_propagator_dt = expm(Hs * -param.dt * 1i);
Ws = eye(2);
Ws(2,2) = -1;
observable = eye(2);
observable(2,2) = -1;
interaction = eye(2);
interaction(2,2) = -1;

% construct BIF-TT
Lb0 = bif_tensor_train(param, method);
Lb = mult_xi_to_Lb(Lb0,param.xi);

% inchoworm solving
FP = inchworm_solve(param, param.M, Lb, interaction, bare_propagator_dt, observable);
rho0 = zeros(2,2);
rho0(1,1) = 1;
idx = @(k, pos) ...
        (k <= -1) * (k + param.N + 1) + ...
        (k >= 1) * (k + param.N + 2) + ...
        (k == 0) * (param.N + 1 + (nargin < 2 || pos));

trace_list = zeros(param.N+1,1);
trace_list(1) = 1;
for n = 1:param.N
    prop = FP{idx(-n),idx(n)};
    ob = trace(rho0 * prop);
    trace_list(n+1) = ob;
end
