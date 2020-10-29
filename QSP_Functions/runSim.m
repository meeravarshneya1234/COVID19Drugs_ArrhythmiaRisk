function datatable = runSim(settings,pert)
%--------------------------------------------------------------------------
                      %% -- runSim.m -- %%
% Description: Runs 2011 O'Hara model, can run either a single cell or
% population

% Inputs:
% --> settings - [struct] simulation protocol settings    
% --> pert - [struct] channel block settings  

% Outputs: 
% --> datatable - [struct] save data for each action potential such as the
% time vector, state variables, and current calculations 
% -------------------------------------------------------------------------
%% 1--- Load Model Parameters 
[p,c] = define_parameters(settings.phys_state,settings.gender);

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
if ~isfield(settings,'scalings')
    if ~isfield(settings,'sigma')
        settings.sigma = 0;
    end
    S = exp(settings.sigma*randn(length(F_G),settings.variations))' ; % same parameter variation for each pop
else
    S = settings.scalings;
end

%% 5--- Define Perturbation Protocol 
[p,c]= perturbations(c,p,pert);
baselineparameters = c;

%% 6--- Run Loop 
parfor ii=1:settings.variations
    scaling = S(ii,:);
    c = scaling_factors(scaling,baselineparameters,F_G);
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


%% Nested Functions 
function [p,c]= perturbations(c,p,pert)
    c.G.GNa  = c.G.GNa  * pert.GNa;
    c.G.GNaL = c.G.GNaL * pert.GNaL;
    c.G.Gto  = c.G.Gto  * pert.Gto;
    c.G.GKr_ = c.G.GKr_ * pert.GKr;
    c.G.GKs_ = c.G.GKs_ * pert.GKs;
    c.G.GK1  = c.G.GK1  * pert.GK1;
    c.G.Gncx = c.G.Gncx * pert.GNCX;
    c.G.GKb  = c.G.GKb  * pert.GKb;
    c.G.GpCa = c.G.GpCa * pert.GpCa;
    c.G.PCa_ = c.G.PCa_ * pert.PCa;
    c.G.Pnak = c.G.Pnak * pert.NaK;
    c.G.PNab = c.G.PNab * pert.PNab;
    c.G.PCab = c.G.PCab * pert.PCab;
    c.G.SERCA_total = c.G.SERCA_total * pert.SERCA ;
    c.G.RyR_total   = c.G.RyR_total   * pert.RyR ;
    c.G.Leak_total  = c.G.Leak_total  * pert.Leak;
    c.G.Trans_total = c.G.Trans_total * pert.Trans;

    % protocols
    p.Ko = p.Ko * pert.Ko;
    p.Inject = pert.Inject;
    p.Cao = p.Cao * pert.Cao;
    p.Nao = p.Nao * pert.Nao;

function c = scaling_factors(scaling,baselineparameters,n_G)
    for iF = 1:length(n_G)
        aF = n_G{iF};
        c.G.(aF) = baselineparameters.G.(aF) * scaling(iF);
    end

