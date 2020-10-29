%% Figure 4D: Virtual population simulations and arrhythmia susceptibility 
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
%% Set Up Simulation Settings 
settings.PCL = 1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 20 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 32; % Stimulus amplitude
settings.steady_state = true; % Start with each scenario's BL model steady state values. 
settings.variations = 1; % run one cell
settings.nBeats = 100; % number of beats to run 
settings.numbertokeep = 2; % determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.gender ='F'; % male - M female -F
settings.phys_state = 'HF'; % healthy or HF
[scalings,labels] = xlsread('COVID-19_Population.csv'); %get population variability values

to_test = [935 395]; % two cells with same GKr,GCaL, but different GNCX
D1 = 7559.6; %concentration of CQ
D2 = 1025.3; %concentration of AZ
Drugs = {'Chloroquine','Azithromycin'};
colors = [0.8500, 0.3250, 0.0980; 0.4940, 0.1840, 0.5560];

for i = 1:length(to_test)
    %% Run PRE-DRUG
    pert = settings_blockcurrents;
    settings.scalings = scalings(to_test(i),:);
    datatable = runSim(settings,pert);
    
    t = datatable.times;
    V = datatable.states(:,1);
    
    figure
    plot(t,V,'k--','linewidth',2)
    hold on        
    
    %% Run DRUG
    pert = settings_blockcurrents;
    pert = get_drug_data(pert,Drugs,[D1;D2]);
    datatable = runSim(settings,pert);
    
    t = datatable.times;
    V = datatable.states(:,1);
    
    plot(t,V,'linewidth',2,'color',colors(i,:))
    hold on
    xlabel('Time (ms)')
    ylabel('V (mV)')
    title(['Cell # ' num2str(to_test(i))])
    set(gca,'FontSize',12,'FontWeight','bold','FontName','Calibri','XGrid','On','YGrid','On')
      
end
    
figure
colormap(colors)
bar(log(scalings(to_test,[4,10,7])'))
xticks(1:3)
xticklabels({'GKr','GCaL','GNCX'})
set(gca,'FontSize',12,'FontWeight','bold','FontName','Calibri','XGrid','On','YGrid','On')

