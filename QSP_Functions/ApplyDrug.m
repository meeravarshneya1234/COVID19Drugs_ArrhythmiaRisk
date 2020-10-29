function ApplyDrug(settings,drugs,Folder)
%--------------------------------------------------------------------------
                            %% -- ApplyDrug.m -- %%
% Description: Apply drug effects on population.  

% Inputs:
% --> settings - [struct] simulation protocol settings 
% --> Folder - [string] main folder to save all the data  
% --> drugs - [struct] drug protocol settings    
%     fieldnames: name - drug names 
%                 fname - folder name to save drug data 
%                 C - drug concentrations 
%--- Example: 
% drug(1).name = {'Lopinavir','Ritonavir'};
% drug(1).fname = 'LP+RT'; % subfolder name 
% drug(1).C = [361.79; 26.35]; % drug concentrations for each drug 

% Outputs: No actual outputs. Data is automatically saved to folder. 
% -------------------------------------------------------------------------
%% Settings
settings.nBeats = 100 ; % Number of beats to simulate
settings.numbertokeep = 10 ;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 0;

%% Run Simulation 
load(fullfile(Folder,'BasePop.mat'),'popICs','popscalings');
[variations,~] = size(popscalings);

n_D = length(drugs);

disp(['Running ' Folder])
for ind = 1:n_D
    D = drugs(ind).name;
    C = drugs(ind).C;
    pert = settings_blockcurrents;
    pert = drugfuncs.calculate_drugblock(pert,D,C);

    % Create folder to save all the data for that trigger
    Filename = drugs(ind).fname;
    disp(['Running ' Filename])
    SubFolder = fullfile(Folder,Filename);
    mkdir(SubFolder);
    
    % Set Up Population Variants
    settings.ICs = popICs;
    settings.scalings = popscalings;
    settings.variations = variations;
    X = runSim(settings,pert); % run simulation

    % Separate cells based on susceptiblity
    Y = popfuncs.find_arrhythmias(settings,X);
    A = sum(Y);
    B = length(Y) - A;
    disp([num2str(A) ' cells formed arrhythmias'])
    
    % Plot 
    figure
    bar(1,B,0.5,'FaceColor',[0, 0.4470, 0.7410])
    hold on
    bar(2,A,0.5,'FaceColor',[0.6350, 0.0780, 0.1840])
    xticks(1:2)
    xticklabels({'(-)arrhythmia','(+)arrhythmia'})
    ylabel('Count','FontSize',14,'FontWeight','bold','FontName','Calibri')
    set(gca,'FontSize',14,'FontWeight','bold','FontName','Calibri','XGrid','on','YGrid','on')
    ylim([0 length(Y)])
    set(gcf,'Position',[520 351 932 747])
    picfile = fullfile(SubFolder, 'TriggerBarPlot.fig');
    saveas(gcf,picfile)
    close(gcf)
    
    % Save Data
    DrugAPs =struct('times',[],'states',[]);
    for i = 1:variations
        [times,states] = popfuncs.splitdata(X(i).times,X(i).states,settings);
        DrugAPs(i).times = times{end};
        DrugAPs(i).states = states{end};
    end
    matfile = fullfile(SubFolder, 'DrugPop.mat');
    save(matfile, 'DrugAPs')
    
    matfile = fullfile(SubFolder, 'Y.mat');
    Y = logical(Y);
    save(matfile, 'Y')
    
    % Clear and reset data for next perturbation
    clear X1 X Y
end
