function rho = ttm_solve(rho0, dynmaps, transfer_tensors, solve_steps)
observable = zeros(1,4); observable(1) = 1; observable(4) = -1;
Kmax = size(dynmaps,2);
rho = cell(1,Kmax + solve_steps + 1);
rho{1} = rho0;
for n = 1: Kmax
    state = zeros(4,1);
    for k = 0 : n-1
        tmp = transfer_tensors{n-k} * rho{k + 1};
        state = state + tmp;
    end
    rho{n+1} = state;
    % rho{n+1} = dynmaps{n} * rho0;
end

for n = Kmax + 1 : Kmax + solve_steps 
    state = zeros(4,1);
    for m = 1 : Kmax
        tmp = transfer_tensors{m} * rho{n - m + 1};
        state = state + tmp;
    end
    rho{n+1} = state;
end


end
