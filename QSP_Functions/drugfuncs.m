classdef drugfuncs
    
    methods(Static)
%--------------------------------------------------------------------------
                            %% -- get_drug_data.m -- %%
% Description: Extract Cmax, IC50, and hill coefficients for the 4 drugs
% based on the Crumb et. al study. See citations. 

% Inputs:
% --> drug - [string] name of drug 
%     options: Azitromycin, Chloroquine, Lopinavir, Ritonavir  

% Outputs:
% --> Cmax - [double] drug Cmax 
% --> IC50 - [double] drug IC50 for 7 channels 
% --> h - [double] drug hill coefficients for 7 channels 
% -------------------------------------------------------------------------
    function [Cmax,IC50,h] = get_drug_data(drug)        
        if strcmp(drug,'Azithromycin')
            Cmax = 1937;
            IC50(1) = NaN;     h(1) = 1;   %ICaL 
            IC50(2) = 70796;   h(2) = 0.5; %IKr
            IC50(3) = NaN;     h(3) = 1;   %IK1
            IC50(4) = 88764;   h(4) = 0.5; %Ito
            IC50(5) = 470131;  h(5) = 1.4; %IKs
            IC50(6) = NaN;     h(6) = 1;   %INa
            IC50(7) = 189128;  h(7) = 1.9; %INaL
        elseif strcmp(drug,'Chloroquine')
            Cmax = 249.5;
            IC50(1) = NaN;     h(1) = 1;
            IC50(2) = 6889;    h(2) = 0.6;
            IC50(3) = 10595;   h(3) = 0.8;
            IC50(4) = NaN;     h(4) = 1;
            IC50(5) = NaN;     h(5) = 1;
            IC50(6) = NaN;     h(6) = 1;
            IC50(7) = NaN;     h(7) = 1;
        elseif strcmp(drug,'Lopinavir')
            Cmax = 703.7;
            IC50(1) = 15601;   h(1) = 1;
            IC50(2) = 5170;    h(2) = 1.2;
            IC50(3) = NaN;     h(3) = 1;
            IC50(4) = NaN;     h(4) = 1;
            IC50(5) = NaN;     h(5) = 1;
            IC50(6) = NaN;     h(6) = 1;
            IC50(7) = NaN;     h(7) = 1;
        elseif strcmp(drug,'Ritonavir')
            Cmax = 436.9;
            IC50(1) = 8228;    h(1) = 1.3;
            IC50(2) = 5157;    h(2) = 1;
            IC50(3) = NaN;     h(3) = 1;
            IC50(4) = NaN;     h(4) = 1;
            IC50(5) = NaN;     h(5) = 1;
            IC50(6) = NaN;     h(6) = 1;
            IC50(7) = 7175;    h(7) = 0.7;
        else
            disp('Choose one of the following: Azithromycin, Chloroquine, Lopinavir, or Ritonavir')
            return 
        end

    end
%--------------------------------------------------------------------------
                      %% -- calculate_drugblock.m -- %%
% Description: Calculate amount of drug block on 7 cardiac channels.  

% Inputs:
% --> pert - [struct] initial channel settings     
% --> drug - [string] drug name; options: Azitromycin, Chloroquine, Lopinavir, Ritonavir  
% --> C - [double] drug concentration 

% Outputs:
% --> pert - [struct] drug block applied to each channel  
% -------------------------------------------------------------------------    
    function pert = calculate_drugblock(pert,drug,C)      
        for i = 1:length(drug) % loop through each drug 
            [~,IC50,h] = drugfuncs.get_drug_data(drug{i}); % extract data 
            DrugBlock(i,:)= IC50.^h./(IC50.^h+C(i).^h); % calculate drug block 
            DrugBlock(i,isnan(DrugBlock(i,:))) = 1; 
        end
        BlockData = prod(DrugBlock,1);
        pert.PCa = pert.PCa*BlockData(1);
        pert.GKr = pert.GKr*BlockData(2);
        pert.GK1  = pert.GK1*BlockData(3);
        pert.Gto = pert.Gto*BlockData(4);
        pert.GKs = pert.GKs*BlockData(5);
        pert.GNa = pert.GNa*BlockData(6);
        pert.GNaL = pert.GNaL*BlockData(7);
        
    end
    
    end 
end

