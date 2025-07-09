function G_list = matrice_to_TT(USV_list,s_len, i, f)
% This function converts a given TPC matrix B to TT form
% Inputs:
%   USV_list = {U,S,V}, the reduced SVD decomposition of B
%   s_len: length of the time sequence s, i.e. M+1
%   (i,f): index for the pair (si,sf)
% Output:
%   G_list: a cell storing TT-cores.

U = USV_list{1}; S = USV_list{2}; V = USV_list{3};
N = size(U,1);
R = size(U,2);
G_list = {};
pvt = 1;
while pvt < i
    if pvt == 1
        % G = zeros(N,R);
        % G(:,1) = ones(N,1);
        G = ones(N,1);
    else
        % G = zeros(R,N,R);
        % G(1,:,1) = ones(1,N);
        G = ones(1,N,1);
    end
    G_list{pvt} = G;
    pvt = pvt + 1;
end %pvt = i
if pvt == 1 %i==1
    G = U;
    G_list{pvt} = U;
else
    % G = zeros(R,N,R);
    % G(1,:,:) = U;
    G = U;
    G = reshape(G,[1,N,R]);
    G_list{pvt} = G;
end
pvt = pvt + 1;%pvt = i+1;
while pvt < f
    G = zeros(R,N,R);
    for s = 1 : N
        for l = 1: R
            for r = 1 : R
                if l == r
                    G(l,s,r) = 1;
                end
            end
        end
    end
    G_list{pvt} = G;
    pvt = pvt + 1;
end %pvt = f
if pvt == s_len
    G = S * V';
    G_list{pvt} = G;
    return;
else
    % G = zeros(R, N, R);
    % G(:,:,1) = S * V';
    G = zeros(R, N, 1);
    G(:,:,1) = S * V';
    G_list{pvt} = G;
end
pvt = pvt + 1; % pvt = f+1;
while pvt < s_len
    % G = zeros(R,N,R);
    % G(1,:,1) = ones(1,N);
    G = ones(1,N,1);
    G_list{pvt} = G;
    pvt = pvt + 1;
end % pvt = s_len
% G = zeros(R,N);
% G(1,:) = ones(1,N);
G = ones(1,N);
G_list{pvt} = G;
end



% 
% function G_list = matrice_to_TT(USV_list,s_len, i, f)
% % This function converts a given TPC matrix B to TT form
% % Inputs:
% %   USV_list = {U,S,V}, the reduced SVD decomposition of B
% %   s_len: length of the time sequence s, i.e. M+1
% %   (i,f): index for the pair (si,sf)
% % Output:
% %   G_list: a cell storing TT-cores.
% 
% U = USV_list{1}; S = USV_list{2}; V = USV_list{3};
% N = size(U,1);
% R = size(U,2);
% G_list = {};
% pvt = 1;
% for n = 1: i-1
%     if n == 1
%         G = zeros(N,R);
%         G(:,1) = ones(N,1);
%     else
%         G = zeros(R,N,R);
%         G(1,:,1) = ones(1,N);
%     end
%     G_list{pvt} = G;
%     pvt = pvt + 1;
% end % pvt = i
% if pvt == 1 % i == 1
%     G = U;
%     G_list{pvt} = U;
% else
%     G = zeros(R,N,R);
%     G(1,:,:) = U;
%     G_list{pvt} = G;
% end
% pvt = pvt + 1; % pvt = i + 1
% 
% for pvt = i+1 : f-1
%     G = zeros(R,N,R);
%     for s = 1 : N
%         for l = 1: R
%             for r = 1 : R
%                 if l == r
%                     G(l,s,r) = 1;
%                 end
%             end
%         end
%     end
%     G_list{pvt} = G;
%     pvt = pvt + 1;
% end % pvt = f
% 
% if pvt == s_len % f == s_len
%     G = S * V';
%     G_list{pvt} = G;
%     return;
% else
%     G = zeros(R, N, R);
%     G(:,:,1) = S * V';
%     G_list{pvt} = G;
% end
% pvt = pvt + 1; % pvt = f + 1;
% while pvt < s_len
%     G = zeros(R,N,R);
%     G(1,:,1) = ones(1,N);
%     G_list{pvt} = G;
%     pvt = pvt + 1;
% end
% G = zeros(R,N);
% G(1,:) = ones(1,N);
% G_list{pvt} = G;
% end


