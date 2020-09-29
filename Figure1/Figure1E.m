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
settings.numbertokeep = 10;
settings.gender ='M';
settings.phys_state = 'healthy';
settings.nBeats = 100;

Cmax_LP = get_drug_Cmax('Lopinavir');
Cmax_RT = get_drug_Cmax('Ritonavir');
hill = 1;

D = logspace(log10(0.3),log10(10),10);
Drugs = {'Lopinavir','Ritonavir'};
for ii = 1:length(D) %for each LP  
    for i = 1:length(D) %for each RT
        pert = settings_blockcurrents;
        CC = [D(ii)*Cmax_LP; D(i)*Cmax_RT];
        pert = get_drug_data(pert,Drugs,CC,hill);
        scalings(i,:) = cell2mat(struct2cell(pert));
        datatable(i) = runSim(settings,pert);
        [time,V,APDs,Cai] = split_data(datatable(i).times,datatable(i).states(:,1),settings,datatable(i).states(:,6));
        APD(ii,i) = find_APD90(time{end},V{end});
    end   
end

% Run Baseline model 
pert_BL = settings_blockcurrents;
settings_BL = settings;
settings_BL.numbertokeep =  1;
datatable_BL = runSim(settings_BL,pert_BL);
APDbase = find_APD90(datatable_BL.times,datatable_BL.states(:,1));

% Get Change in APD 
H = APD-APDbase;
H = flipud(H);

% Create heatmap of change in APD for different drug combinations 
cmap = xlsread('rainbow.xlsx');
figure
imagesc(H)
colormap(cmap)
xticks(1:10)
xticklabels(round(D(1:10),2))
xtickangle(0)
yticks(1:10)
yticklabels(round(fliplr(D),2))
ytickangle(90)
xlabel('RT (nM)')
ylabel('LP (nM)')
set(gca,'FontSize',12,'FontName','Calibri')
colorbar('location','northoutside','box','off','TickDirection','out')
set(gca,'TickDir','out')
set(gcf,'Position',[360 561 343 361])
caxis([10 200])