% tanhFE = 8 * ((Problem.FE-Fes) / (Problem.maxFE-Fes)) - 4;
% K = 0.5-tanh(tanhFE)/2;
% VAR0 = K*mean(cons(index));



% Define the range for Fes
Fes = linspace(1, 100000, 1000);
maxFes = 100000;

% Calculate the values of a
a = 8 * (Fes / maxFes) - 4;

% K=(1-tanh(8 * (Fes / maxFes) - 4))/2;
K=(1-tanh(4* (2*Fes-maxFes)/ maxFes))/2;
% Plot the graph
plot(Fes, K);
xlabel('Iteration');
ylabel('K');

