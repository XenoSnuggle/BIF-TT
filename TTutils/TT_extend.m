function TT_out = TT_extend(TT_in,indices)
%% This function extend the length of TT by length(indices)
% through creating sites at indices
N = size(TT_in{1},1);
len = length(TT_in);
TT_out = cell(1,len+length(indices));
TT_in_idx = 0;
for k = 1:len+length(indices)
    if ismember(k, indices)
        if TT_in_idx == 0
            TT_out{k} = ones(N,1);
        elseif TT_in_idx == len
            TT_out{k} = ones(1,N,1);
        else
            alpha = size(TT_in{TT_in_idx+1},1);
            TT_out{k} = zeros(alpha,N,alpha);
            for n = 1:N
                for j = 1:alpha
                    TT_out{k}(j,n,j) = 1;
                end
            end
        end
    else
        TT_in_idx = TT_in_idx + 1;
        if TT_in_idx == 1
            size_core = size(TT_in{TT_in_idx});
            TT_out{k} = reshape(TT_in{TT_in_idx},[1,size_core]);
        else
            TT_out{k} = TT_in{TT_in_idx};
        end
    end
end
end

