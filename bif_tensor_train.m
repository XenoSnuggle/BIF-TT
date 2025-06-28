function Lb = bif_tensor_train(param, method)
load Extra_Inchworm_Set.mat Extra_Inchworm_Set;
N = param.N;
target_ranks = get_target_rank(param);
disp(["Target ranks = [ ", num2str(target_ranks)," ]."]);
switch method
    case 'tol'
        disp(["Start Building Lb_TT: Tol-Rounding with tol = ", num2str(param.ttrounding_tol)]);
    case 'rank'
        disp(["Start Building Lb_TT: Rank-Rounding with RAM_RANGE = [", num2str(param.ram_range),']']);
    otherwise
        error('Wrong Method!');
end

[Si, Sf] = meshgrid(-N:N, -N:N);
B = tpc_function(param.dt * (abs(Sf) - abs(Si)), param);
% B = param.xi * B;

USV_list = reduced_svd(B, param.svd_tol);

Lb = cell(1,(param.M+1) / 2);
for m = 1:2:param.M
    Lb_m = cell(1,m+1);
    if m == 1
        Lb_m{1} = USV_list{1};
        Lb_m{2} = USV_list{2} * USV_list{3}';
        Lb{1} = Lb_m;
    else
        for k = 3:m
            TT_old = TT_extend(Lb{(m-1)/2},[1,k]);
            TT_new = matrice_to_TT(USV_list,m+1,1,k);
            if k == 3
                Lb_m = TTmult_ver2(TT_old,TT_new);
                clear TT_old TT_new;
            else
                l_old = get_bond_link(TT_old);
                l_new = get_bond_link(TT_new);
                l_after = mult_tt_bond(l_old, l_new);
                if ram_from_bond(l_after) > ((m<=5)*(param.ram_range(1)) + (m>=7) * (param.ram_range(2)))
                    l = bond_divided(l_old,l_new,m);
                    TT_old = tt_rounding(TT_old, l, true);
                end
                tmp_m = TTmult_ver2(TT_old,TT_new);
                clear TT_old TT_new;

                Lb_m = TTadd(Lb_m,tmp_m);
                clear tmp_m;
                if ram_from_bond(get_bond_link(Lb_m)) > ((m<=5)*(param.ram_range(1)) + (m>=7) * (param.ram_range(2)))
                    l = min(get_bond_link(Lb_m), target_ranks((m+1)/2));
                    Lb_m = tt_rounding(Lb_m, l, true);
                end

            end
        end % for k = 3:m
        current_extra = Extra_Inchworm_Set{(m+1)/2};
        for k = 1:size(current_extra,1)
            for l = 1:2:m
                if l == 1
                    G_merge = matrice_to_TT(USV_list,m+1,current_extra(k,l),current_extra(k,l+1));
                else
                    tmp_m = matrice_to_TT(USV_list,m+1,current_extra(k,l),current_extra(k,l+1));
                    l_old = get_bond_link(G_merge);
                    l_new = get_bond_link(tmp_m);
                    l_after = mult_tt_bond(l_old, l_new);
                    if ram_from_bond(l_after) > ((m<=5)*(param.ram_range(1)) + (m>=7) * (param.ram_range(2)))
                        l = bond_divided(l_old, l_new, m);
                        G_merge = tt_rounding(G_merge, l, true);
                    end
                    G_merge = TTmult_ver2(G_merge, tmp_m);
                    clear tmp_m;
                end

            end
            fprintf("Merging Extra (after l-loop, TTadd) for m = %d, k = %d:\n",m,k);
            Lb_m = TTadd(Lb_m,G_merge);
            clear G_merge;
            if ram_from_bond(get_bond_link(Lb_m)) > ((m<=5)*(param.ram_range(1)) + (m>=7) * (param.ram_range(2)))
                l = min(get_bond_link(Lb_m), target_ranks((m+1)/2));
                Lb_m = tt_rounding(Lb_m, l, true);
            end
        end

        Lb{(m+1)/2} = Lb_m;
        fprintf('Finished building Lb_m (after Rounded)m = %d:\n', m);
        fprintf('Lb_m: ');
        print_bond_dim(Lb_m);
        disp(['Ram: ',num2str(get_ram(m,Lb_m)),' GB']);
    end
end



    function target_ranks = get_target_rank(param)
        target_ranks = [];
        disp(['[RAM_MIN, RAM_MAX] = [',num2str(param.ram_range(1)),', ', num2str(param.ram_range(2)),']']);
        get_ram_ = @(n,r,m) ((2*n+1)*r*2 + (2*n+1)*r^2*(m-1))*16/1024^3;
        target_ranks(1) = 0;
        for m = 3:2:param.M
            r = 10000;
            while get_ram_(param.N, r, m) > param.ram_range(1)
                r = r - 100;
            end
            target_ranks((m+1)/2) = r;
            disp(['M = ', num2str(m), ': target rank = ', num2str(r)]);
        end
    end

    function l = get_bond_link(tt)
        d = size(tt,2);
        l = zeros(1,d - 1);
        l(1) = size(tt{1},2);
        for i = 2 : d-1
            l(i) = size(tt{i},3);
        end
    end

    function l = mult_tt_bond(l1,l2)
        l = zeros(1,length(l1));
        for i = 1 : length(l1)
            l(i) =l1(i) * l2(i);
        end
    end

    function ram = ram_from_bond(l)
        N_ = 2* N + 1;
        ram = N_ * l(1);
        for i = 2:length(l)
            ram = ram + N_ * l(i-1) * l(i);
        end
        ram = ram + N_ * l(end);
        ram = 16 * ram/(1024)^3;
    end

    function l = bond_divided(l1, l2, m)
        l = l1;

        mid = ceil(length(l)/2);
        idx = mid;
        for i = 1:mid - 1
            if mid - i >=1
                idx = [idx, mid - i];
            end
            if mid + i <= length(l)
                idx = [idx, mid + i];
            end
        end
        ram = ram_from_bond(mult_tt_bond(l,l2));
        while ram > ((m<=5)*(param.ram_range(1)) + (m>=7) * (param.ram_range(2)))
            if all(l == 10)
                return;
            end
            for i = idx
                if ram < ((m<=5)*(param.ram_range(1)) + (m>=7) * (param.ram_range(2)))
                    return;
                end
                if l(i)>20
                    l(i) = l(i) - 10;
                    ram = ram_from_bond(mult_tt_bond(l,l2));

                end
            end
        end
    end

    function ram = get_ram(m,tt)
        N_ = 2 * N + 1;
        ram = N_ * size(tt{1},2);
        for i = 2:m
            ram = ram + N_ * size(tt{i},1) * size(tt{i},3);
        end
        ram = ram + N_ * size(tt{m+1},1);
        ram = 16 * ram/(1024)^3;
    end

end

