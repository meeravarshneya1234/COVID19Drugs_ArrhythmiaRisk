%% Figure 1A: Concentration-dependent effects of COVID-19 drugs on 
%% ventricular action potentials.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%"Investigational treatments for COVID-19 may increase ventricular
% arrhythmia risk through drug interactions" 

% By:Varshneya,Irurzun-Arana,Campana,Dariolli,Gutierrez,Pullinger,Sobie

% For questions, please contact Dr.Eric A Sobie -> eric.sobie@mssm.edu 
% or put in a pull request or open an issue on the github repository:
% https://github.com/meeravarshneya1234/COVID19Drugs_ArrhythmiaRisk.git. 

%--- Note:
% Results displayed in manuscript were run using MATLAB 2019b on a 64bit
% Intel Processor. For exact replication of figures it is best to use these
% settings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cmax_CQ = get_drug_Cmax('Chloroquine');
Cmax_AZ = get_drug_Cmax('Azithromycin');
Cmax_LP = get_drug_Cmax('Lopinavir');
Cmax_RT = get_drug_Cmax('Ritonavir');
hill = 1; % Use hill coef in drug block equation 
Drugs = {'Azithromycin','Chloroquine','Ritonavir','Lopinavir'};
C = [Cmax_AZ Cmax_CQ Cmax_RT Cmax_LP]; % Cmax values for each drug 
D = 10; % Multiplicative factor for Cmax value

% obtain channel drug block 
for ii = 1:length(Drugs)    
    pert = settings_blockcurrents;
    CC = D*C(ii);
    pert = get_drug_data(pert,Drugs{ii},CC,hill);
    scalings(ii,:) = cell2mat(struct2cell(pert));
end

% create heatmap with drug block 
cmap = readtable('Reds.xlsx','ReadVariableNames', false);
cmap = cmap{:,:};
cmap(1,:) = [1; 1; 1];
figure
imagesc(1-scalings(:,[1:6,10])*100)
xticklabels({'INa','INaL','Ito','IKr','IKs','IK1','ICaL'})
xtickangle(90)
yticks(1:length(Drugs))
yticklabels(Drugs)
ytickangle(0)
set(gca,'FontSize',16,'FontWeight','bold','FontName','Calibri')
colormap(cmap)
colorbar
    