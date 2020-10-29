classdef popfuncs
    
    methods(Static)
%--------------------------------------------------------------------------
                      %% -- create_scale_vector.m -- %%
% Description: creates population variability scaling matrix.  

% Inputs:
% --> settings - [struct] simulation protocol settings    
% --> variations - [double array] number of members of the population 

% Outputs: 
% --> scaling_matrix - [double array] population variability scaling matrix
% -------------------------------------------------------------------------
        function scaling_matrix = create_scale_vector(settings,variations)                     
            [~,c] = define_parameters('healthy','M'); %use this default setting just to get number of parameters to vary
            nG = length(fieldnames(c.G));  
            scaling_matrix = exp(settings.sigma*randn(nG,variations))' ; % same parameter variation for each pop
        end 
               
%--------------------------------------------------------------------------
                      %% -- clean_population.m -- %%
% Description: Determines if there are arrhythmias in the population,
% removes them, and reruns new cells 

% Inputs:
% --> settings - [struct] simulation protocol settings    
% --> Xin - [struct] population time and state variables  

% Outputs: 
% --> Xout - [struct] cleaned population 
% -------------------------------------------------------------------------
        function Xout = clean_population(Xin,settings)
            settings_rerun = settings;
            removes = popfuncs.find_arrhythmias(settings_rerun,Xin);
            X2 = Xin;
            X2(find(removes)) = [];
            n_X2 = length(X2);
            total_variations = settings_rerun.variations;
            
            while n_X2 < total_variations
                settings_rerun.variations = total_variations - n_X2; % number of APs to rerun
                settings_rerun.scalings = popfuncs.create_scale_vector(settings_rerun,settings_rerun.variations);
                disp(['Need to Rerun: ' num2str(settings_rerun.variations) ' cells.'])
                
                pert = settings_blockcurrents;
                X = runSim(settings_rerun,pert); % run population simulation
                X2 = [X2,X];        
                removes = popfuncs.find_arrhythmias(settings_rerun,X2);
                X2(find(removes)) = [];
                n_X2 = length(X2);
            end            
            Xout = X2;
        end
        
%--------------------------------------------------------------------------
                      %% -- find_arrhythmias.m -- %%
% Description: Determines if there are arrhythmias in a population

% Inputs:
% --> settings - [struct] simulation protocol settings    
% --> Xin - [struct] population time and state variables  

% Outputs: 
% --> remove_AP - [logical] logical array of which cells have arrhythmias  
% -------------------------------------------------------------------------
        function remove_AP = find_arrhythmias(settings,Xin)
            remove_AP = false(length(Xin),1);
            for ii = 1:length(Xin) %for each member of the population 
                Xt = Xin(ii).times;
                XV = Xin(ii).states(:,1);
                [times,volts] = popfuncs.splitdata(Xt,XV,settings); % saved the last 10 beats, separate each beat

                % check each beat for arrhythmic activity 
                for i = 1:settings.numbertokeep %for each beat 
                    t = times{i}; t = t - t(1);
                    V = volts{i}(:,1);
                    APDs(i) = find_APD90(t,V);
                    
                    if isnan(APDs(i)) || V(1) > -70 || max(V) < 15 %check if cell failed to repolarize
                        remove_AP(ii) = true;
                         
                    else %check AP for early afterdepolarizations 
                        [~,index] = find(t==settings.stim_delay);
                        time_roi = t(index:end);
                        Vm_roi = V(index:end);
                        dVm_roi = (Vm_roi(2:end)-Vm_roi(1:end-1))./(time_roi(2:end)-time_roi(1:end-1));
                        dVm_dVm = dVm_roi(1:end-1).*dVm_roi(2:end);
                        [~, idx_dVm] = max(dVm_roi);
                        dVm_0 = dVm_dVm(idx_dVm:end) < -1.0*10^-6;
                        dVm_t = time_roi(3:end);
                        tVm_0 = dVm_t(idx_dVm:end);
                        ts = tVm_0(dVm_0);
                        time_between_depols = diff(ts);
                        
                        if any(time_between_depols > 130)
                            remove_AP(ii) = true;  
                        end  
                    end  
                end
                if diff(APDs) > 5 %check for alternans  
                    remove_AP(ii) = true;
                end
            end
%             I = find(remove_AP);
%             Xout(I) = [];

        end
%--------------------------------------------------------------------------
                      %% -- splitdata.m -- %%
% Description: When multiple beats of a single cell are simulated, this
% function separates each beat into its own cell array. This is used
% mainly when settings.numbertokeep is greater than 1. 

% Inputs:
% --> settings - [struct] simulation protocol settings    
% --> Ti - [double array] time matrix 
% --> V - [double array] voltage matrix 

% Outputs: 
% --> times - [cell] time vector for each beat 
% --> volts - [cell] voltage vector for each beat 
% -------------------------------------------------------------------------
        function [times,volts] = splitdata(Ti,V,settings)   
            numberofAPs = settings.numbertokeep;
            PCL = settings.PCL;
            intervals = find(~mod(Ti,PCL));
            times = {};
            volts ={};
            for i = 1:numberofAPs
                times{end+1} = Ti(intervals(i):intervals(i+1));
                volts{end+1} = V(intervals(i):intervals(i+1),:);
            end
        end
    end 
end
