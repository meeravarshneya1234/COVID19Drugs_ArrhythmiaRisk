function [pert] = get_drug_data(pert,drug,C,hill)

if ~exist('hill','var') 
    hill = 0;
end 

T = readtable('Crumb_IC50.xlsx','ReadRowNames',true,'PreserveVariableNames',true);
IC50_data = T{drug,6:12};

if hill == 1
    h_data = T{drug,13:end};
    Norm_IC50= IC50_data.^h_data./(IC50_data.^h_data+C.^h_data);
else
    Norm_IC50= IC50_data./(IC50_data+C);
end
Norm_IC50(isnan(Norm_IC50)) = 1;
BlockData = prod(Norm_IC50,1);

pert.PCa = pert.PCa*BlockData(1);
pert.GKr = pert.GKr*BlockData(2);
pert.GK1  = pert.GK1*BlockData(3);
pert.Gto = pert.Gto*BlockData(4);
pert.GKs = pert.GKs*BlockData(5);
pert.GNa = pert.GNa*BlockData(6);
pert.GNaL = pert.GNaL*BlockData(7);
