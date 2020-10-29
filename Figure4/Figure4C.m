%% Figure 4C: Virtual population simulations and arrhythmia susceptibility 
%% in female heart failure group.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%"Investigational treatments for COVID-19 may increase ventricular
% arrhythmia risk through drug interactions" 

% By: Varshneya & Irurzun-Arana et al. 

% For questions, please contact Dr.Eric A Sobie -> eric.sobie@mssm.edu 
% or put in a pull request or open an issue on the github repository:
% https://github.com/meeravarshneya1234/COVID19Drugs_ArrhythmiaRisk.git. 

%--- Note:
% Results displayed in manuscript were run using MATLAB 2019b on a 64bit
% Intel Processor. For exact replication of figures it is best to use these
% settings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load parameters 
[~,c] = define_parameters('HF','Female');
parameters = table2array(cell2table(struct2cell(c.G)));
[scalings,labels] = xlsread('COVID-19_Population.csv');
allparameters = scalings.*parameters';
X_LOG = log(allparameters);
load(fullfile('HFFemale','Y.mat'));

%% Normalize Data 
[N_trials, N_pars] = size(X_LOG);
for ii=1:N_pars % z-score
    X_LOGISTIC(:,ii)=(X_LOG(:,ii)-mean(X_LOG(:,ii)))/std(X_LOG(:,ii));
end
Y_LOGISTIC = 1-(Y-1); % positive integer!

%% Run Logisitc Regression & Plot B-Matrix 
[B_LOGISTIC,~,stats] = mnrfit(X_LOGISTIC,Y_LOGISTIC);
B = B_LOGISTIC(2:end);

figure
for i = 1:length(B)
    b = bar(i,abs(B(i)));
    if B(i) > 0
        set(b,'Facecolor',[0.6350, 0.0780, 0.1840])
    else
        set(b,'Facecolor',[0, 0.4470, 0.7410])
    end
    hold on
end
title('Probability EAD development')
xticks(1:length(B))
xticklabels(labels(1:length(B)))
xtickangle(45)
set(gca,'XLim',[0 length(B)+1])
set(gca,'FontSize',12,'FontWeight','bold','FontName','Calibri','XGrid','on','YGrid','on')
set(gcf,'Position',[442 534 716 331])