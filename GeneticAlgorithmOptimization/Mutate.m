function [ pNew ] = Mutate(p,mutRate)
% Randomly mutates p at rate mutRate

l = length(p);
pNew = zeros(1,l);
for i = 1:l
    x = rand;
    if x < mutRate
        pNew(i) = 1-p(i);
    else
        pNew(i) = p(i);
    end
end

end

