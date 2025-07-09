function print_bond_dim(Lb)
m = size(Lb,2) - 1;
str = sprintf('bond = [ %d', size(Lb{1},2));
for i = 2:m-1
    str = [str ',' sprintf('%d', size(Lb{i},3))]; 
    % str =[str ','];
end
str = [str ',' sprintf('%d', size(Lb{m},3)) ']'];fprintf('%s\n', str);
