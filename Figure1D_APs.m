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
hill = 1;
figure

% Simulate BL
pert = settings_blockcurrents;
datatable = runSim(settings,pert);
plot(datatable.times,datatable.states(:,1),'k','linewidth',2)
hold on
APD(1) = find_APD90(datatable.times,datatable.states(:,1));

% Simulate 3Cmax AZ
pert = settings_blockcurrents;
C = Cmax_AZ*3;
pert = get_drug_data(pert,{'Azithromycin'},C,hill);
datatable = runSim(settings,pert);
plot(datatable.times,datatable.states(:,1),'linewidth',2)
hold on
APD(2) = find_APD90(datatable.times,datatable.states(:,1));

% Simulate 3Cmax CQ
pert = settings_blockcurrents;
C = Cmax_CQ*3;
pert = get_drug_data(pert,{'Chloroquine'},C,hill);
datatable = runSim(settings,pert);
plot(datatable.times,datatable.states(:,1),'linewidth',2)
hold on
APD(3) = find_APD90(datatable.times,datatable.states(:,1));

% Simulate 3Cmax CQ + 3Cmax AZ
pert = settings_blockcurrents;
C = [Cmax_CQ*3; Cmax_AZ*3];
pert = get_drug_data(pert,{'Chloroquine','Azithromycin'},C,hill);
datatable = runSim(settings,pert);
plot(datatable.times,datatable.states(:,1),'linewidth',2)
hold on
APD(4) = find_APD90(datatable.times,datatable.states(:,1));

xlim([0 500])

legend('bl','AZ','CQ','combined')

fname = 'Figure1D_APs';
folder = 'C:\Users\meera\OneDrive - icahn.mssm.edu\Mount Sinai\Sobie Lab\Thesis\Aim 2\COVID-19\COVID Manuscript\CPT PSP\UpdatedFigures\';
saveas(gcf,[folder fname],'emf')
saveas(gcf,[folder fname],'fig')