%% Figure 1D: Concentration-dependent effects of COVID-19 drugs on 
%% ventricular action potentials.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%"Investigational treatments for COVID-19 may increase ventricular
% arrhythmia risk through drug interactions" 

% By:Varshneya & Irurzun-Arana et al.

% For questions, please contact Dr.Eric A Sobie -> eric.sobie@mssm.edu 
% or put in a pull request or open an issue on the github repository:
% https://github.com/meeravarshneya1234/COVID19Drugs_ArrhythmiaRisk.git. 

%--- Note:
% Results displayed in manuscript were run using MATLAB 2019b on a 64bit
% Intel Processor. For exact replication of figures it is best to use these
% settings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings
%-- Initial Settings 
settings.PCL = 1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 20 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 32; % Stimulus amplitude
settings.steady_state = true; % Start with each scenario's BL model steady state values. 
settings.variations = 1; % run one cell
settings.nBeats = 100; % number of beats to run 
settings.numbertokeep = 1; % determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.gender ='M'; % male - M female -F
settings.phys_state = 'healthy'; % healthy or HF

Cmax_CQ = drugfuncs.get_drug_data('Chloroquine'); % get CQ Cmax 
Cmax_AZ = drugfuncs.get_drug_data('Azithromycin'); % get AZ Cmax 

%% Run Simulation 
%% -- Simulate Baseline
figure
pert = settings_blockcurrents;
datatable = runSim(settings,pert);
plot(datatable.times,datatable.states(:,1),'k','linewidth',2)
hold on
APD(1) = find_APD90(datatable.times,datatable.states(:,1));

%% -- Simulate 3Cmax AZ
Test_Concentration = Cmax_AZ*3; % multiple of Cmax
pert = settings_blockcurrents; % initialize channel block settings
pert = drugfuncs.calculate_drugblock(pert,{'Azithromycin'},Test_Concentration); % calculate drug channel block
datatable = runSim(settings,pert);% run simulation with drug block
plot(datatable.times,datatable.states(:,1),'linewidth',2)
APD(2) = find_APD90(datatable.times,datatable.states(:,1));

%% -- Simulate 3Cmax CQ
Test_Concentration = Cmax_CQ*3; % multiple of Cmax
pert = settings_blockcurrents; % initialize channel block settings
pert = drugfuncs.calculate_drugblock(pert,{'Chloroquine'},Test_Concentration); % calculate drug channel block
datatable = runSim(settings,pert);% run simulation with drug block
plot(datatable.times,datatable.states(:,1),'linewidth',2)
APD(3) = find_APD90(datatable.times,datatable.states(:,1));

%% -- Simulate 3Cmax CQ + 3Cmax AZ
Test_Concentration = [Cmax_CQ*3; Cmax_AZ*3]; % multiple of Cmax
pert = settings_blockcurrents; % initialize channel block settings
pert = drugfuncs.calculate_drugblock(pert,{'Chloroquine','Azithromycin'},Test_Concentration); % calculate drug channel block
datatable = runSim(settings,pert);% run simulation with drug block
plot(datatable.times,datatable.states(:,1),'linewidth',2)
APD(4) = find_APD90(datatable.times,datatable.states(:,1));

xlim([0 500])
legend('Baseline','3*AZ','3*CQ','3*CQ+3*AZ')
xlabel('time (ms)')
ylabel('voltage (mV)')
set(gca,'FontSize',12,'FontWeight','bold','FontName','Calibri','XGrid','On','YGrid','On')