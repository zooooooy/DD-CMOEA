function [Population,Fitness] = EnvironmentalSelection(Population,N,isOrigin)
%Population：当前种群的对象集合。
%N：选择的解决方案数量。
%isOrigin：一个布尔值，指示是否考虑约束条件。

    %% Calculate the fitness of each solution
    if isOrigin 
        Fitness = CalFitness(Population.objs,Population.cons);
    else
        Fitness = CalFitness(Population.objs);
    end

    %% Environmental selection
    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del  = Truncation(Population(Next).objs,sum(Next)-N);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    % Population for next generation
    Population = Population(Next);
    Fitness    = Fitness(Next);
    % Sort the population
    [Fitness,rank] = sort(Fitness);
    Population = Population(rank);
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end
