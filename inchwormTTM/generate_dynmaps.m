function dynmaps = generate_dynmaps(param,Lb,interaction, bare_propagator_dt, threads)

observable = repmat({zeros(2)}, 1, 4);
observable{1}(1,1) = 1;observable{2}(1,2) = 1;observable{3}(2,1) = 1;observable{4}(2,2) = 1;

fp = cell(1,4);

p=parpool(threads);
parfor i = 1:4
    fp{i} = inchworm_solve_dynmap(param, param.M, Lb, interaction, bare_propagator_dt, observable{i});
end
delete(p);

idx = @(k, pos) ...
    (k <= -1) * (k + param.N + 1) + ...
    (k >= 1) * (k + param.N + 2) + ...
    (k == 0) * (param.N + 1 + (nargin < 2 || pos));

rho0 = observable;
dynmaps = cell(1,param.N);
for n = 1:param.N
    dyn = zeros(4);
    for i = 1:4
        for j = 1: 4
            % dyn(i,j) = trace(rho0{j} * fp{i}{idx(-n),idx(n)});
            dyn(j,i) = trace(rho0{j} * fp{i}{idx(-n),idx(n)});
        end
    end
    dynmaps{n} = dyn;
end
end
