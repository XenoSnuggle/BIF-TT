function tt2 = tt_rounding(tt,r,if_RQ)
d = length(tt);
n = size(tt{1},1);

%% RQ decomposition
if if_RQ
    tt1 = cell(1,d);
    [Q,R] = qr(transpose(tt{d}),0);
    tt1{d} = transpose(Q);

    for k = d-1:-1:2
        tmp_core = tns_mult(tt{k},3,R,2);
        tmp_mat  = reshape(tmp_core,[size(tmp_core,1),size(tmp_core,2)*size(tmp_core,3)]);
        [Q,R]    = qr(transpose(tmp_mat),0);
        tt1{k}   = reshape(transpose(Q),[size(Q,2),n,size(Q,1)/n]);
    end

    tt1{1} = tns_mult(tt{1},2,R,2);
else
    tt1 = tt;
end

%% t-SVD
tt2 = cell(1,d);
[U,S,V] = svds(tt1{1},r(1)); 
tt2{1} = U;
for k = 2:d-1
    tmp_core = tns_mult(S*V',2,tt1{k},1);
    tmp_mat  = reshape(tmp_core,[size(tmp_core,1)*size(tmp_core,2),size(tmp_core,3)]);
    [U,S,V]  = svds(tmp_mat,r(k));
    % tt2{k}   = reshape(U,[size(U,1)/n,n,r(k)]);

    if r(k)~=size(U,2) && k == d-1
        tt2{k}   = reshape(U,[size(U,1)/n, n, size(U,2)]);
    else
        tt2{k}   = reshape(U,[size(U,1)/n, n, min(r(k),size(U,2))]);
    end
end

tt2{d} = tns_mult(S*V',2,tt1{d},1);