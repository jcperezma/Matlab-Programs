% clear
% clc
% x = -10:.1:10;
% myFunct = @linearFunction;
% p= [ 4 10 ];
% y = feval(myFunct,x, p,c);
% %y = 4*x+10;
% 
% data(:,1) = x;
% data(:,2) = y;

function [p] = OptimizeGA(myFunct,x,p,c,y)
% y is the data that we are trying to match the model to
% p is the list of parameter values that are found to be the best fit

% This function uses a Genetic Algorithm to determine the best fit of the
% conversion curve

% The optimal solution will be chosen to be the most fit organims in the
% final generation

% Start with creating an initial, random population
n = 3; % number of parameters
numBits = 16; % number of bits 
popSize = 50; % number of organisms created in each generation
pMin = [0.0002 0.0008 0.1]; % parameter minimum values
pMax = [0.25 5 2]; % parameter max values
maxGen = 100; % Max number of iterations
tol = .001; % Margin of error willing to accept
mutRate = .02; % percentage of genes that will be mutated
% End of input
data(:,1) = x;
data(:,2) = y;
% population uses binary vectors to describe values of parameter values
% random starting population
population{1} = randi([0,1],popSize, n*numBits);
gen = 1; % first generation

% Breeding/updating of populations

% Exits if max number of generation has been reached
% Exits if tolerance has been met (lowest prefitness value within
% tolerance specified above)
% Exits if average preFitness does not decrease significantly 4 generations 
% in a row
aveCount = 0; 
% 1 if tolerance met, 0 otherwise
tolMet = 0;
              
while gen <= maxGen && aveCount < 4 && tolMet == 0;
    fprintf('gen %d of %d\n',gen,maxGen)
    % count of offspring created
    oCount = 0;
    % fitness calculation
    % returns new population{gen} that has been reordered from most to
    % least fit
    [preFit{gen}, preFitAve(gen), rawFit{gen}, avRawFit(gen), population{gen}] = ...
        fitnessCalc(population{gen}, pMin, pMax, data,myFunct,c);
    
    % Checking to see if most fit individual is within tolerance
    if preFit{gen}(1) < tol
        tolMet =1;
    end
    
    % Determining number of offspring for each organism of previous generation
    offspring1{gen} = rawFit{gen}/avRawFit(gen);
    offspring2{gen} = floor(offspring1{gen});
    offspring3{gen} = offspring1{gen} - offspring2{gen};
    
    % Creation of next generation
    gen = gen+1;
    % Add in most fit member of previous generation to poulation of 
    population{gen}(1,:) = population{gen-1}(1,:);
    oCount = 1; % incrementing offspring count
    
    % Breeding fit individuals
    i = 1; % individual that is being bred
    while oCount < popSize && i < popSize
        for j = 1:offspring2{gen-1}(i)
            % Breeds with next most fit individuals as determined by the
            % matrix offspring 2
            if 0 < i+j < popSize && 0 < i
                population{gen}(oCount+1,:) = Crossover(population{gen-1}(i,:),...
                    population{gen-1}(i+j,:));
                % keeps track of how many time individuals have already been
                % bred
                offspring2{gen-1}(i+j) = offspring2{gen-1}(i+j)-1;
                % incrementing count to reflect number of current individuals
                % in next generation
                oCount = oCount+1;
            end
        end
        i = i+1;
    end
    
    % If this is not sufficient to populate next generation, the population
    % is re-sorted by which individuals have largest number after decimal
    % for number of offsprings
    if oCount < popSize
        % Sorting most to least fit leftovers
        [offspring3{gen-1}, idx] = sort(offspring3{gen-1});
        % Rearranging population in this order
        pop = population{gen-1}(idx,:);
        % breeds 1 and 2, then 3 and 4 and so on, until the next generation
        % is full populated
        j = 1;
        while oCount < popSize && j < popSize
            population{gen}(oCount+1,:) = Crossover(pop(j,:), pop(j+1,:));
            oCount = oCount +1;
            j = j + 2;
        end
    end
    
    % If next generation is still not fully populated, the last individuals
    % are populated with the best individuals of the previous generation
    population{gen}(oCount+1:popSize,:) = population{gen-1}(oCount+1:popSize,:);
    
    % Next generation is now fully populated
    
    % Mutating next generation, minus best individual from previous
    for i = 2:popSize
        population{gen}(i,:) = Mutate(population{gen}(i,:),mutRate);
    end
end

for i = 1:gen
    pop = population{i}(1,:);
    for j = 1:n
        % Extracting binary vectors for each variable
        binVec = pop((j-1)*numBits+1:j*numBits);
        % calculation of parameter value
        p(1,j) = bin2int(binVec,pMin(j), pMax(j));
    end
    parameters(i,:) = p;
end
parameters;
