function [p,c] = define_parameters(phys_state,gender)
%--------------------------------------------------------------------------
                        %% -- define_parameters.m -- %%
% Description: Extract O'Hara model parameters. Only the endocardial cell
% parameters are used.

% Inputs:
% --> phys_state - [string] physiological state studied; options: healthy or HF 
% --> gender - [string] gender studied; options: 'M' or 'F' 

% Outputs: 
% --> p - [struct] main model parameters 
% --> c - [struct] model parameters that are varied when creating population
% -------------------------------------------------------------------------
%% Standard Parameters 
p.Inject = 0;
p.Nao=140.0;
p.Cao=1.8;
p.Ko=5.4;

%physical constants
p.R=8314.0;
p.T=310.0;
p.F=96485.0;
p.Cm=1.0;   %uF                  

%cell geometry
p.L=0.01;
p.rad=0.0011;
p.vcell=1000*3.14*p.rad*p.rad*p.L;
p.Ageo=2*3.14*p.rad*p.rad+2*3.14*p.rad*p.L;
p.Acap=2*p.Ageo;
p.vmyo=0.68*p.vcell;
p.vnsr=0.0552*p.vcell;
p.vjsr=0.0048*p.vcell;
p.vss=0.02*p.vcell;

%jsr constants
p.bt=4.75;
p.a_rel=0.5*p.bt;

%% Parameters that change based on gender 
if strcmp(gender,'M') ==1
    % computed quantities that do not change during simulation
    c.G.GNa=75;
    c.G.GNaL=0.0075;
    c.G.Gto=0.02;
    c.G.GKr_=0.046;
    c.G.GKs_=0.0034;
    c.G.GK1=0.1908;
    c.G.Gncx=0.0008;
    c.G.GKb=0.003;
    c.G.GpCa=0.0005;
    c.G.PCa_=0.0001;
    c.G.Pnak=30;
    p.cmdnmax=0.05;
    c.G.PNab=3.75e-10;
    c.G.PCab=2.5e-8;
    
    c.G.SERCA_total = 1 ;
    c.G.RyR_total = 1 ;
    c.G.Leak_total = 1;
    c.G.Trans_total = 1;
    p.K_CaMKa = 1;
else 
    c.G.GNa=75;
    c.G.GNaL=0.0075;
    c.G.Gto=0.02;
    c.G.GKr_=0.046;
    c.G.GKs_=0.0034;
    c.G.GK1=0.1908;
    c.G.Gncx=0.0008;
    c.G.GKb=0.003;
    c.G.GpCa=0.0005;
    c.G.PCa_=0.0001;
    c.G.Pnak=30;
    p.cmdnmax=0.05;
    c.G.PNab=3.75e-10;
    c.G.PCab=2.5e-8;   
    c.G.SERCA_total = 1 ;
    c.G.RyR_total = 1 ;
    c.G.Leak_total = 1;
    c.G.Trans_total = 1;
    p.K_CaMKa = 1;
       
    c.G.Gto=c.G.Gto*0.64;
    c.G.GKr_= c.G.GKr_*0.79;
    c.G.GKs_= c.G.GKs_ * 0.83;
    c.G.GK1 = c.G.GK1 *0.86;
    c.G.GpCa= c.G.GpCa *1.6;
    c.G.Pnak= c.G.Pnak * 0.79;
    c.G.SERCA_total = c.G.SERCA_total *1.15;  
end

%% Physiological State: Heart Failure 
if strcmp(phys_state,'HF') ==1

    c.G.GNaL = c.G.GNaL*1.8;
    c.G.Gto = c.G.Gto*0.4;
    c.G.GK1 = c.G.GK1*0.68;
    c.G.Gncx = c.G.Gncx*1.75;
    c.G.Pnak=c.G.Pnak*0.7;
    c.G.SERCA_total = c.G.SERCA_total*0.5 ;
    c.G.Leak_total = c.G.Leak_total*1.3;
    p.K_CaMKa = p.K_CaMKa *1.5;
    
end


