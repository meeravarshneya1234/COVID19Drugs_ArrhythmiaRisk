function datatable = runSim(settings,pert)
%% 1--- Load Model Parameters 
[p,c] = define_parameters(settings.celltype,settings.phys_state,settings.gender);

%% 2--- Load Model Initial Conditions of State Variables 
if ~isfield(settings,'ICs')
    y0 = ICs(settings);
    y0s = repmat(y0,settings.variations,1);
else 
    y0s = settings.ICs;
end

%% 3--- Define Simulation Protocol 
stim_starts = settings.stim_delay + settings.PCL*(0:settings.nBeats-1)  ;
stim_ends = stim_starts + settings.stim_dur ;

% Create intervals for each beat 
simints = 3*settings.nBeats ;
for i=1:settings.nBeats
    intervals(3*i-2,:) = [settings.PCL*(i-1),stim_starts(i)] ; %beginning 
    intervals(3*i-1,:) = [stim_starts(i),stim_ends(i)] ; %stimulus 
    intervals(3*i,:) = [stim_ends(i),settings.PCL*i] ; % stimulus ends 
end
tend = settings.nBeats*settings.PCL ;              % end of simulation, ms
intervals(end,:) = [stim_ends(end),tend] ;

% Determine when to apply stim_amp or 0 amp  
Istim = zeros(simints,1) ;
stimindices = 3*(1:settings.nBeats) - 1 ; % apply stimulus on second part of intervals
Istim(stimindices) = -settings.stim_amp ; 

%% 4--- Population Variables 
F_G = fieldnames(c.G);
F_p = fieldnames(c.p);
F_V = fieldnames(c.V);

if ~isfield(settings,'scalings')
    S1 = exp(settings.sigmaG*randn(length(F_G),settings.variations))' ; % same parameter variation for each pop
    S2 = exp(settings.sigmap*randn(length(F_p),settings.variations))' ; % same parameter variation for each pop
    S3 = (settings.sigmaV*randn(length(F_V),settings.variations))' ; % same parameter variation for each pop
    S = [S1 S2 S3];
else
    S = settings.scalings;
end

%% 5--- Define Perturbation Protocol 
[p,c]= perturbations(c,p,pert);
baselineparameters = c;

%% 6--- Run Loop 
parfor ii=1:settings.variations
    scaling = S(ii,:);
    c = scaling_factors(scaling,baselineparameters,F_G,F_p,F_V);
    statevar_i = y0s(ii,:);
    
    % stimulate cell
    if (settings.nBeats > settings.numbertokeep)
        for i=1:simints-3*settings.numbertokeep
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            [post,posstatevars] = ode15s(@dydt_Ohara,intervals(i,:),statevar_i,options,Istim(i),p,c) ;
            statevar_i = posstatevars(end,:) ;
            t = post(end) ;
        end % for
        statevars = statevar_i ;
        for i=simints-3*settings.numbertokeep+1:simints
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            [post,posstatevars] = ode15s(@dydt_Ohara,intervals(i,:),statevar_i,options,Istim(i),p,c) ;
            t = [t;post(2:end)] ;
            statevars = [statevars;posstatevars(2:end,:)] ;
            statevar_i = posstatevars(end,:) ;
        end % for
    else
        t = 0 ;
        statevars = statevar_i ;
        for i=1:simints
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            [post,posstatevars] = ode15s(@dydt_Ohara,intervals(i,:),statevar_i,options,Istim(i),p,c) ;
            t = [t;post(2:end)] ;
            statevars = [statevars;posstatevars(2:end,:)] ;
            statevar_i = posstatevars(end,:) ;
        end % for
    end % if
    
    outputcell = num2cell(statevars,1) ;
    
    datatable(ii).times =  t - min(t) ;
    datatable(ii).states = statevars;
    datatable(ii).scalings = scaling;
    datatable(ii).currents = calculate_currents(p,c,outputcell);
    datatable(ii).intervals = intervals;


end