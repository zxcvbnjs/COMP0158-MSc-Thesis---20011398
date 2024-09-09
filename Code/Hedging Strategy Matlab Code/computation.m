% Define the path to the Excel file
%filename = '/Users/shiliurui/Documents/MATLAB/Hedging Strategy/simulation_data_2_100000.xlsx';
filename = '/Users/shiliurui/Documents/MATLAB/Hedging Strategy/simulation_data_experiment2_graph_based_approach_100000_points.xlsx';

% Read parameter settings as a matrix since there are no headers
paraSet = readmatrix(filename, 'Sheet', 'ParaSet0'); % Reading directly into a numeric matrix

% Read each sheet into MATLAB as matrices
alphaVec = readmatrix(filename, 'Sheet', 'AlphaVecN'); 
sampleZ = readmatrix(filename, 'Sheet', 'Z_T_Matrix'); 
sampleA = readmatrix(filename, 'Sheet', 'A_T_Matrix');
sampleD = readmatrix(filename, 'Sheet', 'D_T_Matrix');
K=0;
N=length(sampleD);
paraSet(12) = 100000;

%%
[nvCV, nvMVec, nvQVec, selectionVec, nvKappaVec, nvAlphaVec] = getNVEffFront(paraSet, paraSet, alphaVec, sampleD, K, N);
[optMVec,optQVec, optCVVec, optKappaVec, optSelectionVec, optLambdaVec] = getOptEffFront(paraSet,paraSet, alphaVec, sampleA, sampleZ, sampleD, K, N);
%% Draw Payoff Comparison Plot
figure; % Create a new figure
plot(nvAlphaVec, nvCV, 'b'); % 'b' specifies a blue line
hold on; % Keep the current plot so that subsequent plots are added to it
plot(nvAlphaVec, optCVVec, 'r'); % 'r' specifies a red line
xlabel('Alpha Values');
ylabel('CVaR Values'); 
title('Payoff Comparison of Classical NV and Optimal CVaR Model'); % Provide a suitable title
legend('NV CV', 'Optimal CV');
grid on;
hold off;
%saveas(gcf, 'Payoff_Comparison_Plot Graph-Based Experiment.svg'); % Save the figure as SVG
%% Draw Quantity Comparison Plot
figure; % Create a new figure
plot(nvAlphaVec, nvQVec, 'b'); % 'b' specifies a blue line
hold on; % Keep the current plot so that subsequent plots are added to it
plot(nvAlphaVec, optQVec, 'r'); % 'r' specifies a red line
xlabel('Alpha Values');
ylabel('Quantity Values'); 
title('Quantity Comparison of Classical NV and Optimal CVaR Model'); % Provide a suitable title
legend('NV CV', 'Optimal CV');
grid on;
hold off;
%saveas(gcf, 'Quantity_Comparison_Plot.svg'); % Save the figure as SVG
%% Draw Percentage Change Plot
figure; % Create another new figure for percentage change
percentageChange = (optMVec - nvCV) ./ nvCV * 100;
plot(nvAlphaVec, percentageChange, 'r'); % 'g' specifies a red line
xlabel('Alpha Values');
ylabel('Percentage Change (%)');
title('Percentage Change between NV and Optimal Model'); % Provide a suitable title
grid on;
%saveas(gcf, 'Percentage_Change_Plot_Mean-Based_Experiment.png'); % Save the figure as PNG
saveas(gcf, 'Percentage_Change_Plot_Graph-Based_Experiment.png'); % Save the figure as PNG
%%
save('Graph Based Experiment workspace variables.mat');
