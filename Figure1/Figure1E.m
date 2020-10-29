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

%-- Run Baseline Model 
pert_BL = settings_blockcurrents;
settings_BL = settings;
settings_BL.numbertokeep =  1;
datatable_BL = runSim(settings_BL,pert_BL);
APDbase = find_APD90(datatable_BL.times,datatable_BL.states(:,1));

%-- Loop through different combinations of LP+RT
Drugs = {'Lopinavir','Ritonavir'};
Cmax_LP = drugfuncs.get_drug_data('Lopinavir');
Cmax_RT = drugfuncs.get_drug_data('Ritonavir');
D = logspace(log10(0.3),log10(10),10); % 
for ii = 1:length(D) %for each CQ  
    parfor i = 1:length(D) %for each AZ      
        Test_Concentration = [D(ii)*Cmax_CQ; D(i)*Cmax_AZ]; % multiples of Cmax
        pert = settings_blockcurrents; % initialize channel block settings
        pert = drugfuncs.calculate_drugblock(pert,Drugs,Test_Concentration); % calculate drug channel block
        datatable = runSim(settings,pert);% run simulation with drug block
        APD(ii,i) = find_APD90(datatable.times,datatable.states(:,1)); % find APD 
    end   
end

%-- Calculate AP Prolongation 
H = flipud(APD-APDbase);

%-- Create heatmap of change in APD for different drug combinations 
cmap = readmatrix('colormaps.xlsx','Sheet','rainbow');
figure
imagesc(H)
colormap(cmap)
xticks(1:10)
xticklabels(round(D(1:10),2))
xtickangle(0)
yticks(1:10)
yticklabels(round(fliplr(D),2))
ytickangle(90)
ytickangle(90)
xlabel('RT (nM)')
ylabel('LP (nM)')
set(gca,'FontSize',12,'FontName','Calibri')
colorbar('location','northoutside','box','off','TickDirection','out')
set(gca,'TickDir','out')
caxis([10 200])