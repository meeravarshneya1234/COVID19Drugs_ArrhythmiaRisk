%% Figure 1A: Concentration-dependent effects of COVID-19 drugs on 
%% ventricular action potentials.
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

%-- obtain channel drug block 
Drugs = {'Azithromycin','Chloroquine','Ritonavir','Lopinavir'};
D = 10; % study cell under 10*Cmax drug block 
for ii = 1:length(Drugs)
    Cmax = drugfuncs.get_drug_data(Drugs(ii));
    Test_Concentration = D*Cmax;
    pert = settings_blockcurrents; % initialize channel block settings
    pert = drugfuncs.calculate_drugblock(pert,Drugs(ii),Test_Concentration); % calculate drug channel block 
    scalings(ii,:) = cell2mat(struct2cell(pert)); % save channel block to matrix 
end

%-- create heatmap showing channel drug block 
cmap = readmatrix('colormaps.xlsx','Sheet','reds');
figure
imagesc((1-scalings(:,[1:6,10]))*100)
xticklabels({'INa','INaL','Ito','IKr','IKs','IK1','ICaL'})
xtickangle(90)
yticks(1:length(Drugs))
yticklabels(Drugs)
ytickangle(0)
set(gca,'FontSize',16,'FontWeight','bold','FontName','Calibri')
colormap(cmap)
colorbar
