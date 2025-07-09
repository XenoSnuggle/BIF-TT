% function f = TTeval(A,idx)
% %Evaluate a tensor for idx vector
% 
% 
% A{1} = squeeze(A{1});
% 
% f = A{1}(idx(1),:);
% 
% for i=2:numel(A)
%    % disp(['Processing TT: ', num2str(i)]);
%    f = f*squeeze(A{i}(:,idx(i),:));
% end
% end


function f = TTeval(A,idx)
%Evaluate a tensor for idx vector


A{1} = squeeze(A{1});

f = A{1}(idx(1),:);

for i=2:numel(A)
   % disp(['Processing TT: ', num2str(i)]);
   if size(A{i},1) == 1
       f = f*squeeze(A{i}(:,idx(i),:)).';
   else
       f = f*squeeze(A{i}(:,idx(i),:));
   end
end
end
