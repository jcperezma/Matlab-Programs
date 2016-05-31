function [ int ] = bin2int(binVec, mini, maxi)
% Converts binary vector to integer, and then converts integer to variable
% value

% maxi and mini are the maximum and minimu values the parameter in question
% can take

l = length(binVec);
b = 0;
for i = 1:l
    b = b+ binVec(i)*2^(l-i);
end

int = (b/2^l)*(maxi-mini)+mini;

end

