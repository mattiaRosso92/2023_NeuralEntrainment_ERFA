% Quick fix to manually reconstruct bad EEG channels, whenever FieldTrip fails
% to do so; at this moment in time, the approach is not yet scalable: do
% NOT loop over all subjects and conditions, as it works with one electrode
% at a time

% Select your missing channel
missing_chan = {'POz'};
idx_missing = find(strcmpi(missing_chan,interpData(1,1).elec.label));

for condi = 3
    % Insert empty label
    if idx_missing ~= 1
    interpData(condi).label = [interpData(condi).label(1:idx_missing-1) ; {''} ; interpData(condi).label(idx_missing:end)];
    else
    interpData(condi).label = [{''} ; interpData(condi).label];                            
    end
    % Fill-in label
    interpData(condi).label(idx_missing) = missing_chan;

    % Work on timeseries
    % Find indexes of neighbours 
    % TODO: soft-code to find neighbours in structure - problem
    % case sensitive to be solved
    [~,idx_neigh] = ismember(neighbours(idx_missing).neighblabel,interpData(condi).elec.label);

    % Insert empty timeseries
    if idx_missing ~= 1
    interpData(condi).trial{1} = [interpData(condi).trial{1}(1:idx_missing-1,:) ; nan(1,length(interpData(condi).trial{1})) ; interpData(condi).trial{1}(idx_missing:end,:)];            
    else
    interpData(condi).trial{1} = [nan(1,length(interpData(condi).trial{1})) ; interpData(condi).trial{1}];
    end
    % Reconstruct missing timeseries
    interpData(condi).trial{1}(idx_missing,:) = mean(interpData(condi).trial{1}(idx_neigh,:),1);



end


clearvars missing_chan idx_missing idx_neigh idx_resort