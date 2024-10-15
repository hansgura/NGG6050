
%% H. Gura
%Later exercise

%%
subjectTag = 'AB';
dataDirectory = '/Users/gurahs/Downloads/LATERdata/data_mgl/F';
expressCutoff = 0.0;

[data, labels] = later_getData(subjectTag, dataDirectory, expressCutoff);

RTs = data{1}; 
clear data;

%%
laterErrFcn = @(fits) -sum(log(normpdf(RTs, fits(1), fits(2))));

initialGuess = [mean(RTs), std(RTs)];

bestFits = fminsearch(laterErrFcn, initialGuess);

muR_opt = bestFits(1);
deltaS_opt = bestFits(2);

%%
reciprocalRTs = 1 ./ RTs;
initialValues = [mean(reciprocalRTs), std(reciprocalRTs)];

lowerBounds = [0.001, 0.001];  
upperBounds = [1000, 1000]; 

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

[bestFits, fval] = fmincon(laterErrFcn, initialValues, [], [], [], [], lowerBounds, upperBounds, [], options);

muR_opt = bestFits(1);
deltaS_opt = bestFits(2);
%%
opts = optimoptions(@fmincon,    ...
   'Algorithm',   'active-set',  ... 
   'MaxIter',     3000,          ... 
   'MaxFunEvals', 3000);             


problem = createOptimProblem('fmincon',    ...
    'objective',   laterErrFcn,     ...
    'x0',          initialValues,   ... 
    'lb',          lowerBounds,     ... 
    'ub',          upperBounds,     ... 
    'options',     opts);              

gs = GlobalSearch;

[fits(ii,:), nllk] = run(gs, problem);

%%

figure;
histogram(empiricalRTs, 'Normalization', 'probability'); % Empirical data
hold on;
histogram(predictedRTs, 'Normalization', 'probability'); % Model-predicted data
xlabel('Reaction Times (RTs)');
ylabel('Probability');
legend('Empirical', 'Model Prediction');
title('Empirical vs Model-Predicted RT Distribution');
hold off;

%%
residuals = empiricalRTs - predictedRTs;
figure;
plot(residuals);
xlabel('Trial');
ylabel('Residual');
title('Residuals of the Fitted Model');

%%
% We can determine whether the model is a reasonable fit by comparing it to
% the raw data and published literature.