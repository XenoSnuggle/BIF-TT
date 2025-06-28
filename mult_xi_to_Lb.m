function Lb = mult_xi_to_Lb(Lb0,xi)
Lb = Lb0;
for i0 = 1 : size(Lb,2)
    Lb{1,i0}{1,1} = xi^(size(Lb{1,i0},2)/2) *  Lb{1,i0}{1,1};
end
end