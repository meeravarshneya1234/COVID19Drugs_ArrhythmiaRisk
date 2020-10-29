%% Figure 4B: Virtual population simulations and arrhythmia susceptibility 
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
%% Set Up Settings to Apply Drug Effects 
%-- Initial Settings 
settings.PCL = 1000 ;  % Pacing, Interval bewteen stimuli,[ms]
settings.stim_delay = 20 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 32; % Stimulus amplitude

%-- Apply drugs 
drug(1).name = {'Chloroquine','Azithromycin'}; 
drug(1).fname = 'CQ+AZ'; % subfolder name 
drug(1).C = [7559.6; 1025.4]; % drug concentrations for each drug 

drug(2).name = {'Lopinavir','Ritonavir'};
drug(2).fname = 'LP+RT'; % subfolder name 
drug(2).C = [361.79; 26.35]; % drug concentrations for each drug 

%% Run Simulation 
%-- Loop through drug effects on different phenotypes 
genders = {'M','F'};
phys_states = {'healthy','HF'};
Folders = {'HealthyMale','HealthyFemale','HFMale','HFFemale'};
n = 0;
for i = 1:length(phys_states) 
    for ii = 1:length(genders)
        n = n+1;
        settings.gender = genders{ii};
        settings.phys_state = phys_states{i};
        ApplyDrug(settings,drug,Folders{n}); % run simulation 
    end 
end 

%% Plot AP Prolongation 
Drugs = {'CQ+AZ','LP+RT'};
for iii = 1:length(Drugs)
    
    for ii = 1:length(Folders)
        
        load(fullfile(Folders{ii},'BasePop.mat'))
        load(fullfile(Folders{ii},Drugs{iii},'DrugPop.mat'))
        
        for i = 1:length(BaseAPs)
            APD_base(i,ii) = find_APD90(BaseAPs(i).times,BaseAPs(i).states(:,1));
            APD_drug(i,ii) = find_APD90(DrugAPs(i).times,DrugAPs(i).states(:,1));
        end
        
        load(fullfile(Folders{ii},Drugs{iii},'Y.mat'))
        APD_drug(Y,ii) = NaN;
    end
    
    figure
    vs = violinplot(APD_drug-APD_base);
    for i = 1:4
        vs(i).ShowData = 0;
        vs(i).ViolinAlpha = 1;
        vs(i).EdgeColor = [0,0,0];
        vs(i).BoxColor = [0,0,0];
    end
    % Create ylabel
    ylabel('\DeltaAPD','FontSize',13);
    title(Drugs{iii})
    
    % Set the remaining axes properties
    set(gca,'FontName','Calibri','FontSize',12,'GridAlpha',0.1,'GridColor',...
        [0 0 0],'XGrid','on','XTick',[1 2 3 4],'XTickLabel',...
        Folders,'YGrid','on');
    set(gcf,'Position',[360 634 560 288])
    
end
