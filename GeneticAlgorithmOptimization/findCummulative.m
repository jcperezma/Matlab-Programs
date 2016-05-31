function [ accum ] = findCummulative( N_i,l_i )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        Wt = sum(N_i.*l_i);
        Wi = N_i.*l_i /Wt;

        accum(1) = Wi(1);
        for i =2 : length(N_i)
            accum(i) = Wi(i) +accum(i-1);
        end
        
end

