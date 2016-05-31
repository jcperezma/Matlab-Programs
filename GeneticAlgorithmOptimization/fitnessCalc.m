function [ preFitness, preFitAve, rawFitness, avgFitness, popOut ] = ...
    fitnessCalc( pop, pMin, pMax, data,myFunct,c)
% population represents the individuals whose fitness is to be assessed
% pMin and pMax represent the ranges of the parameter values
% data represents both the x and y values of the curve you are trying to
% match
% Calculates the raw fitness, pre fitness  average pre fitness and average 
% fitness of a population, as well as returns the population sorted from
% most to least fit

x = data(:,1);
matrix_C = data(:,2);
% Extracting binary vectors for each variable
[h,l] = size(pop); % h = number of individuals
n = length(pMin); % number of variables
numBin = l/n; % size of binomial vectors
preFitness = zeros(h,1); % Matrix to hold pre fitness values
p = zeros(n,1);
for i = 1:h
    for j = 1:n
        % Extracting binary vectors for each variable
        binVec = pop(i,(j-1)*numBin+1:j*numBin);
        % calculation of parameter value
        p(j) = bin2int(binVec,pMin(j), pMax(j));
    end
    % calculation of model
    %[y,yPrime] = ComputeCure(p(1), p(2), p(3), p(4), p(5), data)
    y = feval(myFunct,x', p,c);
    %y = p(1)*x+p(2);
    % prefitness (sum of squares of difference)
    % preFitness(i) = sum((matrix_C(:,2)-y).^2);
    preFitness(i) = sum((matrix_C-y').^2);
end

% Calculation of raw fitness
preFitAve = mean(preFitness);
% Sorting most to least fit
[preFitness, idx] = sort(preFitness);
% Rearranging population in this order
popOut = pop(idx,:);
% least fit individual
maxP = max(preFitness); 
rawFitness = maxP-preFitness; 
avgFitness = mean(rawFitness);

end
