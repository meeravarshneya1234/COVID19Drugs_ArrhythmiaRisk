%% Figure 3C: Sex differences and cardiac disease influence cardiac effects
%% of COVID-19 drugs. 
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
%-- Initial Settings 
settings.PCL = 1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 20 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 32; % Stimulus amplitude
settings.sigma = 0; % no variability
settings.steady_state = true; % start with steady state values 
settings.variations = 1; % run one cell
settings.nBeats = 100; % number of beats to run 
settings.numbertokeep = 1; % determine how many beats to keep. 1 = last beat, 2 = last two beats

%-- Drug Settings 
Drugs = {'Lopinavir','Ritonavir'};
Test_Concentration = [361.79; 26.35];% concentration of LP and RT (from PK results Fig 2)

%-- Loop through four phenotypes to study drug effects on AP
gender ={'M','F','M','F'};
phys_state = {'healthy','healthy','HF','HF'};
figure
for ii = 1:length(gender)
    settings.gender = gender{ii};
    settings.phys_state = phys_state{ii};
    
    % Run simulation WITHOUT drug 
    subplot(1,4,ii)
    pert = settings_blockcurrents;
    datatable = runSim(settings,pert);
    plot(datatable.times,datatable.states(:,1),'k','linewidth',2)
    hold on
    APD(1,ii) = find_APD90(datatable.times,datatable.states(:,1));
    
    % Run simulation WITH drug
    pert = settings_blockcurrents; % initialize channel block settings
    pert = drugfuncs.calculate_drugblock(pert,Drugs,Test_Concentration); % calculate drug channel block
    datatable = runSim(settings,pert);% run simulation with drug block
    plot(datatable.times,datatable.states(:,1),'linewidth',2)
    APD(2,ii) = find_APD90(datatable.times,datatable.states(:,1));
    
    strs = {strcat(gender{ii}, ':' ,phys_state{ii}), strcat('\DeltaAPD=' ,num2str(round(APD(2,ii)- APD(1,ii),2)) ,'ms')};
    title(strs)
    box(gca,'on');
    set(gca,'FontName','Calibri Light','FontSize',12,'XGrid','on',...
        'XMinorTick','on');
    xlim([0 500])
end
set(gcf,'Position',[2291 198.3333 1.8573e+03 420.0000])