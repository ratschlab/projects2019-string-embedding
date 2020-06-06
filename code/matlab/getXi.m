function x = getXi(Data, i)
if iscell(Data)
    x = Data{i};
else
    x = Data(i,:);
end