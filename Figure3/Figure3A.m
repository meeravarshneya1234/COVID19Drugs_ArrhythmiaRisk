%% Figure 3A: Sex differences and cardiac disease influence cardiac effects
%% of COVID-19 drugs. 
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
genders = {'M','F'};
phys_states = {'healthy','HF'};
headers = {};
parameters = [];

% Loop through the different phenotypes 
for i = 1:length(phys_states)    
    for ii = 1:length(genders)
        phys_state = phys_states{i}; % healthy or HF
        gender = genders{ii}; % M vs F
               
        [p,c] = define_parameters(phys_state,gender);
        F_G = fieldnames(c.G);
        parameter_labels = F_G(:)';
        parameter_labels{end+1} = 'CamKa';
                
        parameters(:,end+1) = [cell2mat(struct2cell(c.G));p.K_CaMKa];
        headers{end+1} = strcat(genders{i},':',phys_states{ii});
    end
end

% create matrix of parameters 
T = array2table(parameters);
T.Properties.VariableNames = headers;
T.Properties.RowNames = parameter_labels;
params_norm = parameters./parameters(:,1); %normalize parameters to healthy male 

figure
bar(log10(params_norm([2:7,9,11,14,16,18],:)))
xticks(1:length([2:7,9,11,14,16,18]))
xticklabels(parameter_labels([2:7,9,11,14,16,18]))
xtickangle(45)
ylabel('Parameter Values (Log Scale)')
xlabel('Parameters')
legend(headers)
set(gca,'FontName','Calibri Light','FontSize',12,'XGrid','on','XMinorTick','on');
set(gcf,'Position',[2.4537e+03 76.3333 1.4153e+03 580])