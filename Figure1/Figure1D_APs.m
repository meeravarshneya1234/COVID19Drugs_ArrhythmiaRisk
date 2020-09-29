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
settings.nBeats = 100;

Cmax_CQ = get_drug_Cmax('Chloroquine');
Cmax_AZ = get_drug_Cmax('Azithromycin');
hill = 1;
figure

% Simulate BL
pert = settings_blockcurrents;
datatable = runSim(settings,pert);
plot(datatable.times,datatable.states(:,1),'k','linewidth',2)
hold on
APD(1) = find_APD90(datatable.times,datatable.states(:,1));

% Simulate 3Cmax AZ
pert = settings_blockcurrents;
C = Cmax_AZ*3;
pert = get_drug_data(pert,{'Azithromycin'},C,hill);
datatable = runSim(settings,pert);
plot(datatable.times,datatable.states(:,1),'linewidth',2)
hold on
APD(2) = find_APD90(datatable.times,datatable.states(:,1));

% Simulate 3Cmax CQ
pert = settings_blockcurrents;
C = Cmax_CQ*3;
pert = get_drug_data(pert,{'Chloroquine'},C,hill);
datatable = runSim(settings,pert);
plot(datatable.times,datatable.states(:,1),'linewidth',2)
hold on
APD(3) = find_APD90(datatable.times,datatable.states(:,1));

% Simulate 3Cmax CQ + 3Cmax AZ
pert = settings_blockcurrents;
C = [Cmax_CQ*3; Cmax_AZ*3];
pert = get_drug_data(pert,{'Chloroquine','Azithromycin'},C,hill);
datatable = runSim(settings,pert);
plot(datatable.times,datatable.states(:,1),'linewidth',2)
hold on
APD(4) = find_APD90(datatable.times,datatable.states(:,1));

xlim([0 500])
legend('bl','AZ','CQ','combined')
set(gca,'FontSize',12,'FontWeight','bold','FontName','Calibri','XGrid','On','YGrid','On')

