function drug_Cmax = get_drug_Cmax(drug)

[data,labels] = xlsread('Crumb_IC50.xlsx');
drug_names = labels(2:end,1);
col_names = labels(1,2:end);
d_index = ismember(drug_names,drug);
c_index = strcmp('Free Cmax',col_names);
drug_Cmax = data(d_index,c_index);
