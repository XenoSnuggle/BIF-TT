function fp = inchworm_solve_dynmap(param, M, Lb, interaction, bare_propagator_dt, observable)
NLb = (size(Lb{1,1}{1,1},1) - 1)/2;
N = param.N;
dt = param.dt;
% Hs = zeros(2,2);
% Hs(1,1) = param.epsilon; Hs(2,2) = -param.epsilon;
% Hs(1,2) = param.Delta; Hs(2,1) = param.Delta;
% bare_propagator_dt = expm(Hs * -param.dt * 1i);
% interaction = eye(2);
% interaction(2,2) = -1;
% observable = eye(2);
% observable(2,2) = -1;
identity = eye(2);
fp = cell(2*N+2, 2*N+2);
for i = 1: 2 * N + 2
    fp{i,i} = identity;
end

for i = N - 1: -1: 0
    solve_step(-N, -i, M);
    fp{i+N + 2,2*N + 2} = fp{ 1,-i+N + 1}'; %G(5,6) = G(-6,-5)'
    % obj.solve_step(i, N, M);
    for k = -i + N + 2 : N + 1
        fp{k - (N - i),k} = fp{1,-i + N + 1}; % G(-5,-4)=G(-6,-5)
    end
    for k = i + N + 2 : -1 : N + 2 % N-1 to 0+
        fp{k,k + N - i} = fp{i + N + 2,2*N + 2}; % G(4,5) = G(5,6)
    end
end

fp{N+1, N+2} = observable;

for i = 1:N
    fp{N+1, i+N+2} = fp{N +2, i + N + 2} * observable;
    fp{N+1-i, N+2} = observable * fp{N+1-i, N+1};
end

for i = 1:N
    j = i;
    for j = i:N
        solve_step(-i,j,M); %G(-i,j), i = 1:N, j=i:N
        % fp{-j+N + 1,i+N + 2} = fp{-i+N + 1,j+N + 2}'; %G(-2,1) = G(-1,2)'
        if j~=i
            solve_step(-j,i,M);
        end
    end
end


    function solve_step(mi, mf, M)
        disp(['Solving mi = ', num2str(mi), ', mf = ', num2str(mf)]);
        if mi <0 && mf > 0
            cross = true;
        else
            cross = false;
        end

        if cross
            mi_idx = mi + N + 1;
            mf_idx = mf + N + 2;
        elseif mi >= 0
            mi_idx = mi + N + 2;
            mf_idx = mf + N + 2;
        elseif mf <= 0
            mi_idx = mi + N + 1;
            mf_idx = mf + N + 1;
        end


        G0 = fp{mi_idx, mf_idx-1};
        G_star = G0;

        fp{mi_idx, mf_idx} = G_star;
        for m = 1: 2: M
            integral = integrate(mi, mf-1, cross, m);
            G_star = G_star + dt * integral;
            fp{mi_idx, mf_idx} = G_star;
        end

        if mf > 0
            G_star = conj(bare_propagator_dt) * G_star;
        else
            G_star = bare_propagator_dt * G_star;
        end

        fp{mi_idx, mf_idx} = G_star;

        G1 = G_star;
        for m = 1: 2: M
            integral = integrate(mi, mf, cross, m);
            G1 = G1 + dt * integral;
        end
        if mf > 0
            G = conj(bare_propagator_dt) * G0 + G1;
        else
            G = bare_propagator_dt * G0 + G1;
        end
        fp{mi_idx, mf_idx} = G./2;
    end


    function integral = integrate(mi, mf, cross, m)
        if mi == mf
            integral = 0;
            return
        end

        Lb_m = Lb{(m+1)/2};

        if cross
            n = mf - mi + 1;
            if mf == 0
                n = mf - mi;
            end
        else
            n = mf - mi;
        end
        weight = @(k,i,f) (k==i || k==f) * 0.5 + (k~=i && k~=f) * 1;

        idx = @(k, pos) ...
            (k <= -1) * (k + N + 1) + ...
            (k >= 1) * (k + N + 2) + ...
            (k == 0) * (N + 1 + (nargin < 2 || pos));

        if m==1
            integral = integrate_m2(Lb_m,cross);

        else
            I = integrate_I_first(Lb_m{m},Lb_m{m+1});
            for mm = 2:m-1
                I = integrate_I_mid(I,Lb_m{m-mm+1});
            end
            integral = integrate_I_final(I, Lb_m{1});
        end

        function integral = integrate_m2(Lb_2,cross)
            integral = zeros(2,2);
            if cross
                for k = mi : 0
                    % sgn = -1;
                    wei = weight(k,mi,0);
                    G1 = fp{idx(k,0),idx(mf,1)};
                    G0 = fp{idx(mi,0),idx(k,0)};
                    bif = TTeval(Lb_2,[k + NLb + 1, mf + NLb + 1]);
                    tmp = wei * interaction * G1 * ...
                        interaction * G0 * bif;
                    integral = integral - (1i)^2 * (dt) * tmp;
                end
                if mf >0
                    for k = 0 : mf
                        % sgn = 1;
                        wei = weight(k,0,mf);
                        G1 = fp{idx(k,1),idx(mf,1)};
                        G0 = fp{idx(mi,0),idx(k,1)};
                        bif = TTeval(Lb_2,[k + NLb + 1, mf + NLb + 1]);
                        tmp = wei * interaction * G1 * ...
                            interaction * G0 * bif;
                        integral = integral + (1i)^2 * (dt) * tmp;

                    end
                end

            else
                for k = mi : mf
                    % sgn = 1;
                    wei = weight(k,mi,mf);
                    G1 = fp{idx(k,(mi>=0)),idx(mf,(mi>=0))};
                    G0 = fp{idx(mi,(mi>=0)),idx(k,(mi>=0))};
                    bif = TTeval(Lb_2,[k + NLb + 1, mf + NLb + 1]);
                    tmp = wei * interaction * G1 * ...
                        interaction * G0 * bif;

                    integral = integral + (1i)^2 * (dt) * tmp;
                end
            end

        end

        function I_new = integrate_I_first(L1,L2)

            r1 = size(L1,1);
            r2 = size(L1,3);
            omega_new = zeros(2,2,n+1,r1);
            I_new = zeros(2,2,n,r1);

            if cross
                for alpha1 = 1: r1
                    for k2 = mi:0
                        L = 0;
                        for alpha2 = 1:r2
                            L = L + L1(alpha1,k2 + NLb + 1,alpha2) ...
                                * L2(alpha2, mf + NLb + 1);
                        end
                        G2 = fp{idx(k2,0), idx(mf,1)};
                        omega_new(:,:,k2 - mi + 1, alpha1) = ...
                            interaction * G2 * ...
                            interaction * L;
                    end
                    for k2 = 0: (mf==0)*(-1)+~(mf==0)*mf
                        L = 0;
                        for alpha2 = 1:r2
                            L = L + L1(alpha1,k2 + NLb + 1,alpha2) ...
                                * L2(alpha2, mf + NLb + 1);
                        end
                        G2 = fp{idx(k2,1), idx(mf,1)};
                        omega_new(:,:,k2 - mi + 2, alpha1) = ...
                            interaction * G2 * ...
                            interaction * L;
                    end
                    for k1 = mi : min(0,mf-1)
                        prop = zeros(2,2);
                        for k2 = k1 : 0
                            if k1 == 0
                                break;
                            end
                            wei = weight(k2,k1,0);
                            G1 = fp{idx(k1,0), idx(k2,0)};
                            prop = prop - wei * ...
                                omega_new(:,:,k2 - mi + 1, alpha1) * G1;
                        end
                        for k2 = 0 : (mf==0)*(-1)+~(mf==0)*mf
                            wei = weight(k2,0,mf);
                            G1 = fp{idx(k1,0), idx(k2,1)};
                            prop = prop + wei * ...
                                omega_new(:,:,k2 - mi + 2, alpha1) * G1;
                        end
                        I_new(:,:,k1 - mi + 1, alpha1) = prop;
                    end

                    for k1 = 0: mf - 1
                        prop = zeros(2,2);
                        for k2 = k1 : mf
                            wei = weight(k2, k1, mf);
                            G1 = fp{idx(k1,1), idx(k2,1)};
                            prop = prop + wei * ...
                                omega_new(:,:,k2 - mi + 2, alpha1) * G1;
                        end
                        I_new(:,:,k1 - mi + 2, alpha1) = prop;
                    end
                end

            else
                for alpha1 = 1:r1
                    for k2 = mi : mf
                        L = 0;
                        for alpha2 = 1:r2
                            L = L + L1(alpha1,k2 + NLb + 1,alpha2) ...
                                * L2(alpha2, mf + NLb + 1);
                        end
                        G2 = fp{idx(k2,(mi>=0)), idx(mf,(mi>=0))};
                        omega_new(:,:,k2 - mi + 1, alpha1) = ...
                            interaction * G2 * ...
                            interaction * L;
                    end
                    for k1 = mi : mf - 1
                        prop = zeros(2,2);
                        for k2 = k1 : mf
                            wei = weight(k2, k1, mf);
                            G1 = fp{idx(k1,(mi>=0)), idx(k2,(mi>=0))};
                            prop = prop + wei * ...
                                omega_new(:,:,k2 - mi + 1, alpha1) * G1;
                        end
                        I_new(:,:,k1 - mi + 1, alpha1) = prop;
                    end
                end

            end
        end



        function I_new = integrate_I_mid(I_old, L)
            r1 = size(L,1);
            r2 = size(L,3);
            % I_new = zeros(2,2,n,r1);
            omega_new = zeros(2,2,n,r1);
            I_new = zeros(2,2,n,r1);
            if cross
                for alpha1 = 1:r1
                    for k2 = mi : min(mf-1,0)
                        F_slice = squeeze(I_old(:, :, k2 - mi + 1, :)); %[ 2, 2, r2]
                        L_slice = L(alpha1, k2 + NLb + 1, :);     % [1, 1, r2]
                        omega_new(:,:,k2-mi+1,alpha1) = ...
                            sum(F_slice .* reshape(L_slice, [1, 1, r2]), 3) ...
                            * interaction;
                    end

                    for k2 = 0 : (mf==0)*(-1)+ ~(mf==0)*mf-1
                        F_slice = squeeze(I_old(:, :, k2 - mi + 2, :)); %[ 2, 2, r2]
                        L_slice = L(alpha1, k2 + NLb + 1, :);     % [1, 1, r2]
                        omega_new(:,:,k2 - mi + 2, alpha1) = ...
                            sum(F_slice .* reshape(L_slice, [1, 1, r2]), 3) ...
                            * interaction;
                    end
                    for k1 = mi : min(mf-1,0)
                        prop = zeros(2,2);
                        for k2 = k1: min(mf-1,0)
                            if k1 == 0
                                break;
                            end
                            wei = weight(k2, k1, 0);
                            % sgn = -1;
                            G1 = fp{idx(k1,0),idx(k2,0)};
                            prop = prop - wei * ...
                                omega_new(:,:,k2-mi+1,alpha1) * G1;
                        end
                        for k2 = 0 : (mf==0)*(-1)+ ~(mf==0)*mf-1
                            wei = weight(k2, 0, mf);
                            % sgn = 1;
                            G1 = fp{idx(k1,0),idx(k2,1)};
                            prop = prop + wei * ...
                                omega_new(:,:,k2-mi+2,alpha1) * G1;
                        end
                        I_new(:,:,k1 - mi + 1, alpha1) = prop;
                    end
                    for k1 = 0 : mf - 1
                        prop = zeros(2,2);
                        for k2 = k1 : mf -1
                            wei = weight(k2, k1, mf);
                            % sgn = 1;
                            G1 = fp{idx(k1,1),idx(k2,1)};
                            prop = prop + wei * ...
                                omega_new(:,:,k2-mi+2,alpha1) * G1;
                        end
                        I_new(:,:,k1 - mi + 2, alpha1) = prop;
                    end
                end
            else
                for alpha1 = 1:r1
                    for k2 = mi:mf-1
                        F_slice = squeeze(I_old(:, :, k2 - mi + 1, :)); %[ 2, 2, r2]
                        L_slice = L(alpha1, k2 + NLb + 1, :);     % [1, 1, r2]
                        omega_new(:,:,k2 - mi + 1,alpha1) = ...
                            sum(F_slice .* reshape(L_slice, [1, 1, r2]), 3) ...
                            * interaction;
                    end

                    for k1 = mi : mf - 1
                        prop = zeros(2,2);
                        for k2 = k1 : mf - 1

                            wei = weight(k2, k1, mf);
                            % sgn = 1;

                            G1 = fp{idx(k1,(mi>=0)),idx(k2,(mi>=0))};
                            prop = prop + wei * omega_new(:,:,k2-mi+1,alpha1) * G1;
                        end
                        I_new(:,:,k1-mi + 1,alpha1) = prop;
                    end
                end


            end
        end

        function I = integrate_I_final(I_old,L)
            I = zeros(2,2);
            r = size(L,2);
            omega_new = zeros(2,2,n);

            if cross
                for k = mi : min(mf-1,0)
                    F_slice = squeeze(I_old(:, :, k - mi + 1, :)); %[ 2, 2, r]
                    L_slice = L(k + NLb + 1, :);     % [1, r]
                    omega_new(:,:,k - mi + 1) = ...
                        sum(F_slice .* reshape(L_slice, [1, 1, r]), 3) ...
                        * interaction;

                end
                for k = 0 : (mf==0)*(-1)+ ~(mf==0)*mf-1
                    F_slice = squeeze(I_old(:, :, k - mi + 2, :)); %[ 2, 2, r]
                    L_slice = L(k + NLb + 1, :);     % [1, r]
                    FL = sum(F_slice .* reshape(L_slice, [1, 1, r]), 3);
                    omega_new(:,:,k - mi + 2) = ...
                        sum(F_slice .* reshape(L_slice, [1, 1, r]), 3) ...
                        * interaction;
                end
                for k = mi : min(mf-1,0)
                    wei = weight(k,mi,0);
                    G = fp{idx(mi,0),idx(k,0)};
                    I = I - wei * omega_new(:,:,k-mi+1) * G;
                end

                for k = 0 : (mf==0)*(-1)+ ~(mf==0)*mf-1
                    G = fp{idx(mi,0),idx(k,1)};
                    wei = weight(k,0,mf);
                    I = I + wei * omega_new(:,:,k-mi+2) * G;
                end
            else
                for k = mi : mf - 1
                    F_slice = squeeze(I_old(:, :, k - mi + 1, :)); %[ 2, 2, r]
                    L_slice = L(k + NLb + 1, :);     % [1, r]
                    omega_new(:,:,k - mi + 1) = ...
                        sum(F_slice .* reshape(L_slice, [1, 1, r]), 3) ...
                        * interaction;
                end
                for k = mi : mf - 1
                    if mi >= 0
                        G = fp{idx(mi,1), idx(k,1)};
                    elseif mf <= 0
                        G = fp{idx(mi,0), idx(k,0)};
                    end
                    wei = weight(k,mi,mf);
                    I = I + wei * omega_new(:,:,k-mi+1) * G;
                end

            end
            I = (1i)^(m+1) * (dt)^m * I;
        end


    end

end