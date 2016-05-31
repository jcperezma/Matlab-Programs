function [ p3 ] = Crossover(p1,p2)
% Creates offspring (p3) of two organisms (p1,p2) using crossover
c = 1;
% c denotes type of crossover

l = length(p1);
p3 = zeros(1,l);
if c == 1;
    % double point crossover
    xA = randi([1,l]);
    xB = randi([1,l]);

    if xA > xB
        x1 = xB;
        x2 = xA;
    else
        x1 = xA;
        x2 = xB;
    end

    p3(1:x1) = p1(1:x1);
    p3(x1+1:x2) = p2(x1+1:x2);
    p3(x2+1:l) = p1(x2+1:l);
end

if c == 2 
    % Random crossover
    for i = 1:l
        x = rand;
        if x < .5
            p3(l) = p1(l);
        else
            p3(l) = p2(l);
        end
    end
end

            
            


end

