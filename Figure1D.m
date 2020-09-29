%% Figure 1D: Concentration-dependent effects of COVID-19 drugs on 
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
settings.celltype = 'endo'; % 'epi', 'endo', 'mid',
settings.PCL = 1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 20 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 32; % Stimulus amplitude
settings.sigmaG = 0;
settings.sigmap = 0;
settings.sigmaV = 0;
settings.steady_state = 1; % if set to 0, use the manuscript values or provide with own set of ICs
settings.variations = 1;
settings.numbertokeep = 1;
settings.gender ='M';
settings.phys_state = 'healthy';
settings.nBeats = 1;

Cmax_CQ = get_drug_Cmax('Chloroquine');
Cmax_AZ = get_drug_Cmax('Azithromycin');
hill = 1;

D = logspace(log10(0.3),log10(10),10);
Drugs = {'Chloroquine','Azithromycin'};
for ii = 1:length(D) %for each CQ  
    for i = 1:length(D) %for each AZ
        pert = settings_blockcurrents;
        CC = [D(ii)*Cmax_CQ; D(i)*Cmax_AZ];
        pert = get_drug_data(pert,Drugs,CC,hill);
        datatable = runSim(settings,pert);
        APD(ii,i) = find_APD90(datatable.times,datatable.states(:,1));
    end   
end

% Run Baseline model 
pert_BL = settings_blockcurrents;
settings_BL = settings;
settings_BL.numbertokeep =  1;
datatable_BL = runSim(settings_BL,pert_BL);
APDbase = find_APD90(datatable_BL.times,datatable_BL.states(:,1));

% Get Change in APD 
H = APD-APDbase;
H = flipud(H);

% Create heatmap of change in APD for different drug combinations 
cmap = xlsread('rainbow.xlsx');
figure
imagesc(H)
colormap(cmap)
xticks(1:10)
xticklabels(round(D(1:10),2))
xtickangle(0)
yticks(1:10)
yticklabels(round(fliplr(D),2))
ytickangle(90)
xlabel('AZ (nM)')
ylabel('CQ (nM)')
set(gca,'FontSize',12,'FontName','Calibri')
colorbar('location','northoutside','box','off','TickDirection','out')
set(gca,'TickDir','out')
set(gcf,'Position',[360 561 343 361])
caxis([10 200])