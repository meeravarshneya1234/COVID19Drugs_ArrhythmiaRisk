%% Figure 4A: Virtual population simulations and arrhythmia susceptibility 
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

%% Universal Settings for Building Populations 
%-- Initial Settings 
settings.PCL = 1000 ;  % Pacing, Interval bewteen stimuli,[ms]
settings.stim_delay = 20 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 32; % Stimulus amplitude

%-- Population Settings
settings.sigma = 0.3; % standard deviation to vary conductances 
settings.variations = 200; % population size 
settings.remove_arrhythmias = true; % remove cells that form arrhythmias before drug application and replace them with new ones

% To create the population scaling matrix, uncomment the following line: 
% settings.scaling_matrix = popfuncs.create_scale_vector(settings,settings.total_variations); % create population scaling matrix                    
% After creating the matrix, we saved it to a 'COVID-19_Population.csv' file and called that
% going forward. 
settings.scalings = xlsread('COVID-19_Population.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Healthy Male Population 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Initial Settings 
settings.gender = 'M'; % male - M female -F
settings.phys_state = 'healthy'; % healthy or HF

%-- Folder to Save Data 
Folder = 'HealthyMale';
mkdir(Folder)

%-- Run Population 
pert = settings_blockcurrents; 
BuildPopulation(settings,Folder,pert)
set(gcf,'Name','Healthy Male');% Note: manuscript figure only displays 20 of 1000 for simplicity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Healthy Female Population 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Initial Settings 
settings.gender = 'F'; % % male - M female -F
settings.phys_state = 'healthy'; % healthy or HF

%-- Folder to Save Data 
Folder = 'HealthyFemale';
mkdir(Folder)

%-- Run Population 
pert = settings_blockcurrents; % no perturbations, so use the defaults.
BuildPopulation(settings,Folder,pert)
set(gcf,'Name','Healthy Female');% Note: manuscript figure only displays 20 of 1000 for simplicity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run HF Male Population 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Initial Settings 
settings.gender = 'M'; % male - M female -F
settings.phys_state = 'HF'; % healthy or HF

%-- Folder to Save Data 
Folder = 'HFMale';
mkdir(Folder)

%-- Run Population 
pert = settings_blockcurrents; % no perturbations, so use the defaults.
BuildPopulation(settings,Folder,pert)
set(gcf,'Name','HF Male');% Note: manuscript figure only displays 20 of 1000 for simplicity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run HF Female Population 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Initial Settings 
settings.gender = 'F'; % male - M female -F
settings.phys_state = 'HF'; % healthy or HF

%-- Folder to Save Data 
Folder = 'HFFemale';
mkdir(Folder)

%-- Run Population 
pert = settings_blockcurrents; % no perturbations, so use the defaults.
BuildPopulation(settings,Folder,pert)
set(gcf,'Name','HF Female');% Note: manuscript figure only displays 20 of 1000 for simplicity 