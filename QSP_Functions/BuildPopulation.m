function BuildPopulation(settings,Folder,pert)
%--------------------------------------------------------------------------
                        %% -- BuildPopulation.m -- %%
% Description: Build population

% Inputs:
% --> settings - [struct] simulation protocol settings 
% --> Folder - [string] main folder to save all the data  
% --> pert  - [struct] initialize channel block settings 

% Outputs: No actual outputs. Data is automatically saved to folder. 
% -------------------------------------------------------------------------
%% Settings
settings.nBeats = 500; % Number of beats to simulate
settings.numbertokeep = 10 ;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = true; % Start with each scenario's BL model steady state values. 

%% Run Simulation

% Create a Folder to store all of the population data. 
PopFolder = [Folder '\Population\'];
mkdir(PopFolder)
disp('Running New Population...')

% Does the user provide a parameter matrix or should we make one? 
if ~isfield(settings,'scaling_matrix') % make new one
    settings.scalings = popfuncs.create_scale_vector(settings,settings.variations);
end

% Create Population
X = runSim(settings,pert); % run simulation

% Clean matrix to remove arrhythmias that occur pre-trigger
if settings.remove_arrhythmias 
    Xnew = popfuncs.clean_population(X,settings);
else
    Xnew = X;
end

% Save the ICs and scalings
popICs = [];
popscalings = [];
len = length(Xnew);
for i = 1:len
    popICs(end+1,:) = Xnew(i).states(end,:);
    popscalings(end+1,:) = Xnew(i).scalings;
end

%% Run Population using final steady state ICs 
settings2 = settings;
settings2.nBeats = 1 ; % Only one beat
settings2.numbertokeep = 1 ;% Save only that one beat
settings2.steady_state = false; % don't use previous SS, use the new values saved
settings2.ICs = popICs;
settings2.scalings = popscalings;

% Run each cell in population with new steady state ICs 
BaseAPs = runSim(settings2,pert); % run simulation

figure
for ii = 1:settings2.variations
    plot(BaseAPs(ii).times,BaseAPs(ii).states(:,1),'linewidth',2)
    hold on
end
set(gca,'FontSize',14,'FontWeight','bold','FontName','Calibri','XGrid','on','YGrid','on')
xlabel('time (ms)','FontSize',14,'FontWeight','bold','FontName','Calibri')
ylabel('Voltage (mV)','FontSize',14,'FontWeight','bold','FontName','Calibri')
picfile = fullfile(Folder, 'APs.png');
saveas(gcf,picfile)

matfile = fullfile(Folder, 'BasePop.mat');
save(matfile, 'BaseAPs','popICs','popscalings')

