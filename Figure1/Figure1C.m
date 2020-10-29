%% Figure 1C: Concentration-dependent effects of COVID-19 drugs on 
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

%-- Drug Settings 
D = logspace(log10(0.3),log10(10),20); % Multiplicative factor for Cmax value 
Drugs = {'Lopinavir','Ritonavir','Azithromycin','Chloroquine'};
figure
for i = 1:length(Drugs)
    drug = Drugs{i};
    Cmax = drugfuncs.get_drug_data(drug); % get drug's Cmax 
    parfor ii = 1:length(D)        
        Test_Concentration = D(ii)*Cmax; % multiples of Cmax  
        pert = settings_blockcurrents; % initialize channel block settings
        pert = drugfuncs.calculate_drugblock(pert,drug,Test_Concentration); % calculate drug channel block
        datatable = runSim(settings,pert);% run simulation with drug block 
        
        t = datatable.times;
        V = datatable.states(:,1);
        APD(i,ii) = find_APD90(t,V);
    end
    semilogx(D,APD(i,:),'Marker','.','MarkerSize',20,'LineWidth',2)
    hold on
end
xlabel('Multiples of EFTPC')
ylabel('APD (ms)')
xlim([0.3 10])
set(gca,'XMinorTick','on','XScale','log','XTick',[0.3 0.5 1 2 5 10],...
    'XTickLabel',{'0.3','0.5','1','2','5','10'});
legend(Drugs);
set(gca,'FontSize',12,'FontWeight','bold','FontName','Calibri','XGrid','On','YGrid','On')



