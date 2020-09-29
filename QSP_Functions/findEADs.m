function [EAD_index] = findEADs(X1,settings)

F = fieldnames(X1);
len = length(X1.(F{1}));

EAD_index =zeros(len,1);

for ii = 1:len
    Xt = X1.times{ii};
    XV = X1.states{ii,1}(:,1);
    if settings.numbertokeep > 1
        [time,volts,APDs] = split_data(Xt,XV,settings); % split data based on number of beats kept
    else 
        time = {Xt};
        volts = {XV};
        APDs = find_APD90(Xt,XV);
    end   
        
    [APfails,nEADs] = cleandata(APDs,time,volts,settings);
    [~,ind_failed] = find(APfails ==1); % number of failed to repolarize
    [~,ind_EADs] = find(nEADs==1); % number of EADs
    total = [ind_EADs ind_failed];
       
    if ~isempty(total) %no EADs
        EAD_index(ii) = 1;
    end
    
end
EAD_index = logical(EAD_index);