function transfer_tensors = generate_transfer_tensors(dynmaps)
Kmax = size(dynmaps, 2);
transfer_tensors = cell(1,Kmax);
for n = 1 : Kmax
    T = dynmaps{n};
    for m = 1:n - 1
        tmp = transfer_tensors{n-m} * dynmaps{m};
        T = T - tmp;
    end
    transfer_tensors{n} = T;
end
end