%% Figure 1C: Concentration-dependent effects of COVID-19 drugs on 
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
Cmax_LP = get_drug_Cmax('Lopinavir');
Cmax_RT = get_drug_Cmax('Ritonavir');
hill = 1;
D = logspace(log10(0.3),log10(10),20);
Cmax = [Cmax_LP Cmax_RT Cmax_AZ Cmax_CQ];
Drugs = {'Lopinavir','Ritonavir','Azithromycin','Chloroquine'};

figure
fig = gcf;

for i = 1:length(Drugs)
    for ii = 1:length(D)
        % Drug
        pert = settings_blockcurrents;
        CC = D(ii)*Cmax(i);
        pert = get_drug_data(pert,Drugs(i),CC,hill);
        datatable = runSim(settings,pert);
        
        t = datatable.times;
        V = datatable.states(:,1);
        Cai = datatable.states(:,6);
        [outputs,labels] = calculate_features(V,Cai,t);
        APD(i,ii) = outputs(7);

    end
    figure(fig);
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



