function [C] = TTmult(A,B)
%Multiplying two tensor train (non-periodic)
% C = cell(numel(A),1);
C = cell(1,numel(A));

for i=1:numel(A)

    if i~=1&&i~=numel(A)
        n = size(A{i},2);
        r1 = size(A{i},1);
        r2 = size(A{i},3);
        t1 = size(B{i},1);
        t2= size(B{i},3);
        
        tmp = zeros(r1 * t1,n,r2 * t2);
    elseif i==1
        tmp = zeros(size(A{i},1),size(A{i},2)*size(B{i},2));
        n   = size(A{i},1);
    elseif i==numel(A)
        tmp = zeros(size(A{i},1)*size(B{i},1),size(A{i},2));
        n   = size(A{i},2);
    end

    for k=1:n
        if i~=1&&i~=numel(A)
            Ak = squeeze(A{i}(:,k,:));
            Bk = squeeze(B{i}(:,k,:)); 
            if r1 == 1
                 Ak = Ak.';
            end
            if t1 == 1
                Bk = Bk.';
            end
            tmp(:,k,:) = kron(Ak,Bk);
        elseif i==1

            tmp(k,:) = kron(A{i}(k,:), B{i}(k,:));
        elseif i==numel(A)
            tmp(:,k)   = kron(A{i}(:,k),B{i}(:,k));
        end
    end
    C{i} = tmp;
end


