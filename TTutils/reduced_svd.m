function USV_list = reduced_svd(B, tol)
[U_,S_,V_] = svd(B,"econ");
singular_values = diag(S_);
indices = singular_values > tol;
U = U_(:, indices);
S = S_(indices, indices);
V = V_(:, indices);
USV_list = cell(1,3);
USV_list{1} = U;USV_list{2} = S; USV_list{3} = V;
end