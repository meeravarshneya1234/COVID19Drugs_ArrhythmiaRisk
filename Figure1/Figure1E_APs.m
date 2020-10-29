%% Figure 1E: Concentration-dependent effects of COVID-19 drugs on 
%% ventricular action potentials.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%"Investigational treatments for COVID-19 may increase ventricular
% arrhythmia risk through drug interactions" 

% By:Varshneya and Irurzun-Arana et al.

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
settings.sigma = 0; % no variability
settings.steady_state = true; % Start with each scenario's BL model steady state values. 
settings.variations = 1; % run one cell
settings.nBeats = 100; % number of beats to run 
settings.numbertokeep = 1; % determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.gender ='M'; % male - M female -F
settings.phys_state = 'healthy'; % healthy or HF

%-- Drug Settings 
Cmax_LP = drugfuncs.get_drug_data('Lopinavir');
Cmax_RT = drugfuncs.get_drug_data('Ritonavir');

%% Run Simulation 
%% -- Simulate Baseline
figure
pert = settings_blockcurrents;
datatable = runSim(settings,pert);
plot(datatable.times,datatable.states(:,1),'k','linewidth',2)
hold on
APD(1) = find_APD90(datatable.times,datatable.states(:,1));

%% -- Simulate 3Cmax LP
Test_Concentration = Cmax_LP*3; % multiple of Cmax
pert = settings_blockcurrents; % initialize channel block settings
pert = drugfuncs.calculate_drugblock(pert,{'Lopinavir'},Test_Concentration); % calculate drug channel block
datatable = runSim(settings,pert);% run simulation with drug block
plot(datatable.times,datatable.states(:,1),'linewidth',2)
APD(2) = find_APD90(datatable.times,datatable.states(:,1));

%% -- Simulate 3Cmax RT
Test_Concentration = Cmax_RT*3; % multiple of Cmax
pert = settings_blockcurrents; % initialize channel block settings
pert = drugfuncs.calculate_drugblock(pert,{'Ritonavir'},Test_Concentration); % calculate drug channel block
datatable = runSim(settings,pert);% run simulation with drug block
plot(datatable.times,datatable.states(:,1),'linewidth',2)
APD(3) = find_APD90(datatable.times,datatable.states(:,1));

%% -- Simulate 3Cmax LP + 3Cmax RT
Test_Concentration = [Cmax_LP*3; Cmax_RT*3]; % multiple of Cmax
pert = settings_blockcurrents; % initialize channel block settings
pert = drugfuncs.calculate_drugblock(pert,{'Lopinavir','Ritonavir'},Test_Concentration); % calculate drug channel block
datatable = runSim(settings,pert);% run simulation with drug block
plot(datatable.times,datatable.states(:,1),'linewidth',2)
APD(4) = find_APD90(datatable.times,datatable.states(:,1));

xlim([0 500])
legend('Baseline','3*LP','3*RT','3*LP+3*RT')
xlabel('time (ms)')
ylabel('voltage (mV)')
set(gca,'FontSize',12,'FontWeight','bold','FontName','Calibri','XGrid','On','YGrid','On')