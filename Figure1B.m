%% Figure 1B: Concentration-dependent effects of COVID-19 drugs on 
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
settings.gender = 'M';
settings.phys_state = 'healthy';

% Drug Block Data 
Cmax_CQ = get_drug_Cmax('Chloroquine');
Cmax_AZ = get_drug_Cmax('Azithromycin');
Cmax_LP = get_drug_Cmax('Lopinavir');
Cmax_RT = get_drug_Cmax('Ritonavir');
C = [Cmax_AZ Cmax_CQ Cmax_RT Cmax_LP];
Drugs = {'Azithromycin','Chloroquine','Ritonavir','Lopinavir'};
D = 10;
hill = 1; % use hill coef 

[colors,num,typ] = brewermap(8,'Set1');
map = [colors(4,:);colors(2,:);colors(3,:);colors(5,:);colors(1,:)];

% Run baseline model without drug block 
pert = settings_blockcurrents;
settings.nBeats = 1;
datatable = runSim(settings,pert);
tbase = datatable.times;
Vbase = datatable.states(:,1);
APD_base = find_APD90(tbase,Vbase);

% Loop through different drugs 
for ii = 1:length(Drugs)
    figure
    plot(tbase,Vbase,'k','linewidth',2)
    hold on
    settings.nBeats = 100;
    pert = settings_blockcurrents;
    CC = D*C(ii);
    pert = get_drug_data(pert,Drugs(ii),CC,hill);
    datatable = runSim(settings,pert);
    
    t = datatable.times;
    V = datatable.states(:,1);
    Cai = datatable.states(:,6);
    APD(ii) = find_APD90(t,V);
    
    figure(gcf)
    plot(t,V,'linewidth',2,'color',map(ii,:))
    hold on
    title([Drugs{ii} ': ' num2str(APD(ii) - APD_base)])
    set(gca,'FontSize',12,'FontWeight','bold','FontName','Calibri','XGrid','On','YGrid','On')
    
    figure(gcf)
    xlim([0 500])
    ylim([-100 50])

end
    

