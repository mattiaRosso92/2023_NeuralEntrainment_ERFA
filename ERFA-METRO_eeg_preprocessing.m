%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Event-related frequency adjustment (ERFA).
% Preprocessing pipeline
%
% Up to RESS, we treat all the conditions as the same.
% From RESS onwards, we only analyse Tap conditions: 
% loop over Tap conditions and produce 2 gedData outputs. SAVE
% We will separately analyze ERFA trials in Tap conditions
% 
%
% Mattia Rosso , Ghent 16/4/21
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

% Set Fieldtrip path
%https://www.fieldtriptoolbox.org/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path/
restoredefaultpath
addpath '/Users/mattiaipem/Documents/MATLAB/fieldtrip-20201205'
ft_defaults



%% Import data

% Selection for loops
sub = 14;                        % pick subject
% Set directory
path = ['/Users/mattiaipem/Desktop/ERI_' num2str(sub) '_eeg'];
cd (path);

% Settings
srate = 1000; %sampling rate (1kHz, i.e. 1 datapoint every 1ms)

% Conditions are already sorted (upon export from eego software)
% Note: all 'Listen' conditions come first
condlabels = {'Listen_Metro','Listen_Instr1','Listem_Instr2', 'Listen_Lyrics1','Listen_Lyrics2', ...
    'Tap_Metro_Phase','Tap_Metro_Tempo','Tap_Instr_Phase','Tap_Instr_Tempo','Tap_Lyrics_Phase','Tap_Lyrics_Tempo'}; 
nconds = length(condlabels);
lstncond = [1:floor(nconds/2)];
tapconds = [floor(nconds/2)+1:nconds];
timeconds = [60 60 60 60 60 465 465 465 465 465 465]; %in seconds; 'Listen lasts 1 minute only'
ntrials = 40; % hard code, even if there might be 41 sometimes, due to randomization of perturbation-free window length

% Initialize data structure (concatenate at the end of the loop)
macroData = struct([]); 
cfgTrial = struct([]);

for condi = 7%[1 6 7]  
       
    % Macro epoch (from Max 'start' to 'stop')
    cfg=[]; % refresh configuration 
    cfg.dataset = ['ERI_' num2str(sub) '_Condition' num2str(condi) '.vhdr']; %header to read in
    cfg.trialdef.eventtype = 'Start'; % ask FieldTrip by setting '?'
    cfg.trialdef.eventvalue = 's1'; % 's1' for all perturbations via BNC
    % NOTE: at one point, we have to drop first 3 secs. Best if we do it
    % already. DOUBLE-CHEK what happens in the experiment;
    % UPDATE 17/5/22: we are not dropping first 3 seconds anymore
    cfg.trialdef.prestim  = 0;   % in seconds
    cfg.trialdef.poststim = timeconds(condi); % New default in Max application 
    cfg.trialdef.ntrials  = 1;   % takes the first event ("Start")   
    % Import raw timestamps
    cfg.traw = csvread(['ERI_' num2str(sub) '_Condition' num2str(condi) '.csv']);
    cfg.trialfun = 'ft_trialfun_ERI_study'; % custom instead of 'ft_trialfun_general'    
    % Defines trials in dedicated configuration; needed before ft_preprocessing,
    % but will apply way later in the script (after actusal pre-processing)  
    trialTemp = ft_definetrial(cfg); % store trial info  in temporary structure
    %trialTemp = rmfield(trialTemp,'progress'); % remove superfluous field to match structures
    cfgTrial = [cfgTrial , trialTemp] ; % concatenate 
    
    % Import continuous data (no cut, no alignment to 'time 0')
    cfg = [];
    cfg.continous = 'yes';
    cfg.dataset = ['ERI_' num2str(sub) '_Condition' num2str(condi) '.vhdr']; %header to read
    cfg.channel = 'all';
    cfg.demean = 'yes';   % remove DC offset (by mean subtraction in the frequency domain)
    cfg.implicitref = 'CPz';
    % recover CPz; comment for checking the position
    % of the electrode.   
    % Process continuous data 
    contData = ft_preprocessing(cfg);
    
    
    % Filter data
    cfg = [];
    cfg.bsfilter = 'yes'; %notch
    cfg.bsfreq = [49 51; 99 101; 149 151];
    cfg.hpfilter = 'yes'; %high-pass filter for slow drifts
    cfg.hpfreq = 1;
    cfg.lpfilter = 'yes'; %low-pass filter for hi-frex muscular activity
    cfg.lpfreq = 40;
    filtData = ft_preprocessing(cfg , contData);
    
    % Check that the experimenter did not set by mistake
    % "BIPOLAR" amplifier instead of 'Single'. Would it be the case, 
    % remove here extra BIP chans (66:end) before re-arranging order     
    if size(filtData.trial{1},1) > 65
        bipidx = strfind(filtData.label,'BIP'); %find indexes of labels containing 'BIP'
        filtData.trial{1}(~cellfun(@isempty,bipidx),:) = []; % remove BIP timeseries
        filtData.label(~cellfun(@isempty,bipidx)) = []; % remove BIP labels
    end
       
    % Re-arrange electrodes to match layout
    % Store EOG in separate cell; to further remove from channels
    eog = filtData.trial{1}(32,:);
    %exclude here EOG
    cfg              = [];
    cfg.channel      = {'all' '-EOG'};
    filtData   = ft_selectdata(cfg, filtData);
    %move CPz in its position (49), and Oz in (31) 
    tempLabel = [filtData.label(1:30);filtData.label(end-1);filtData.label(31:47);filtData.label(end);filtData.label(48:end-2)];
    tempData  = [filtData.trial{1}(1:30,:);filtData.trial{1}(end-1,:);filtData.trial{1}(31:47,:);filtData.trial{1}(end,:);filtData.trial{1}(48:end-2,:)];
    filtData.label = tempLabel;
    filtData.trial{1} = tempData;
    clear tempLabel;
    clear tempData;

    % Segment into "macro trial" (i.e., from 'time 0' to end) 
    % (use previous trial configuration and continuous filtered data as
    % input arguments)
    macroTemp = ft_redefinetrial(trialTemp, filtData); % temporary assignment
           
    % Visual inspection 
    cfg.trialdef.prestim = 0;
    cfg.viewmode  = 'vertical';
    cfg.blocksize = 8;
    ft_databrowser(cfg,macroTemp); 
    pause;

    % Concatenare structures: https://www.mathworks.com/hel p/matlab/matlab_prog/concatenate-structures.html
    macroData = [macroData , macroTemp] ; % concatenate [contData tempData] ; check along which dimension   
   
end
clear macroTemp;
clear trialTemp
clear eog;


% Define layout, based on Waveguard electrodes locations
cd ('/Users/mattiaipem/Documents/IPEM/Projects/ Lou _ Entrainment/Event_Related_Instability');
cfg = [];
cfg.elec = ft_read_sens('standard_waveguard64.elc'); % read electrodes from document (3D coordinates)
cfg.rotate = 90; % correct the 90deg rotation
layout = ft_prepare_layout(cfg); % this will be used for ICA visualization only
% Visualize layout
ft_layoutplot(cfg);  


%% Visual rejection bad channels
% (NB: Visual rejection bad trials will be performed on GED component)

%If you have a long recording, you will notice things tend to get worse
%towards the end, due to gel drying and electrodes losing contact; however,
%there is a sweetspot with ANT-Neuro were impedances tend to improve within the first hour of recording 

%ALWAYS CHECK for electrodes with variance 0 which are not the reference!! (i.e., 'flat'). They MUST be
%removed, or the ICA computation is compromised

% Initialize data structure 
cleanData = struct([]); 
badChan = cell(1,nconds);

% Loop over conditions (in every block)
for condi = 1:length(macroData)
    
% Refresh cfg
cfg = [];
cfg.method = 'summary'; %(?)
cfg.alim = 5e-5;     %set limits for the graph
cleanTemp = ft_rejectvisual(cfg,macroData(condi));
pause;

% Log bad channels and trials, note indicator of interest e.g. variance; 
badChan{condi} = macroData(condi).label( ~ismember(macroData(condi).label, cleanTemp.label) );

%Rereference AFTER bad chan removal!!! So the eccessive noise does not
%leak into the average; 
cfg = [];
cfg.reref='yes';
cfg.refchannel={'all'}; 

cleanTemp = ft_preprocessing(cfg , cleanTemp);


% Inspection after re-referencing: expect more homogeneous variance
cfg.viewmode = 'vertical';
ft_databrowser(cfg,cleanTemp);
disp(['You are inspecting subject Condition #' num2str(condi)])
pause;



% Concatenare structures: https://www.mathworks.com/hel p/matlab/matlab_prog/concatenate-structures.html
cleanData = [cleanData cleanTemp];


end

% Clear workspace
clear contData; 
clear filtData;
clear macroData;
clear cleanTemp;


%% Compute ICA

% Data should be as similar as possible; it is optimal to run 1 condition
% at a time.

cfg = [];
% With 'runica'
cfg.method       = 'runica';
cfg.channel      = {'all' '-CPz'}; % remove reference, for being a linear combination of all channels
cfg.trials       = 'all';  
cfg.numcomponent = 'all';
cfg.demean       = 'no';

% Initialize ICA cell
components = cell(nconds,1);
for condi = 1:length(cleanData)
% Note this calculation is not deterministic (if method is not JADE): it might sliiightly 
% vary if you iterate; anyway, it's not so important in practice

components{condi} = ft_componentanalysis(cfg, cleanData(condi));

end


% Save dataset so far
% Set directory
cd (path);
save(['ERI_' num2str(sub) '_ICAmatrix_metronomes'] , '-v7.3');


% Check-point
close all;
clear all;


%% ICA artifact removal (by visual inspection)

clear all;
close all;
clc


% ALWAYS check that you are working with the right set!!
% Selection for loops
sub = 20;                        % pick subject
% Set directory
path = ['/Users/mattiaipem/Desktop']; % move to Desktop; since the last Mac update, I cannot access files on HDD
cd (path);

% Import data with ICA matrix
load(['ERI_' num2str(sub) '_ICAmatrix_metronomes']); % check that the matrix is the right one for the current conditions

% Inspect indipendent components
% Loop over conditions 
for condi = 1%:length(cleanData)
    
    % Plot the components topography and activation timeseries
    cfg = [];
    cfg.component = 1:round(length(components{condi}.label)/2);       % specify the component(s) that should be plotted
    cfg.layout    = layout; % specify the layout file that should be used for plotting (and just for plotting)
    cfg.comment   = 'no';
    ft_topoplotIC(cfg, components{condi})
    % Further inspection of the time course of the components
    cfg = [];
    cfg.layout = layout; % specify the layout file that should be used for plotting
    cfg.viewmode = 'component';
    cfg.blocksize = 8;
    ft_databrowser(cfg, components{condi});

    % Annotate components below during pause
    pause;
    
    
end

%% Remove the bad components and backproject the data

% Insert picked bad components, for every condition
         
badComp = {  [  ]; ...
 [  ]; ...
 [  ]};


% Initialize structure         
postIcaData = struct([]);
% Reject components 
for condi = 1:length(cleanData)
    
    cfg = [];
    cfg.component = badComp{condi}; % to be removed component(s)
    % Assign to temporary variable
    postIcaTemp = ft_rejectcomponent(cfg, components{condi}, cleanData(condi));
    % Concatenate structure
    postIcaData = [postIcaData postIcaTemp];

end
clear postIcaTemp

%% Compare data before and after 

% Consider inspecting a smll cluster of frontal electrodes, 
% with a vector of equally spaced offsetts to add to the timeseries,
% to get the 'vertical' view

cond2inspect = [1 2 3]; % select condition to inspect; make it a loop with pause()
% Select channel and trial 
chan2plot = 'Fz';  % Fpz ideal for assessing blink removal; Fz ideal for periodic pad-related artifact
trialIdx = 1; % select trial to inspect
exInterval = 15000:37000;   % select your time window

% Visualize
figure;
for condi = cond2inspect
    chanIdx = find(strcmp(postIcaData(condi).label, chan2plot));
    
    subplot(3,1,condi)
    plot(cleanData(condi).trial{trialIdx}(chanIdx,exInterval));
    hold on
    plot(postIcaData(condi).trial{trialIdx}(chanIdx,exInterval) , 'r');
    legend('Before ICA removal' , 'After ICA removal');
    ylim([-100 100]);
    title(['Channel ' num2str(chan2plot)  ' - Trial ' num2str(trialIdx) ' - Condition #1']);

end


%% Bad channels reconstruction (Interpolation)

% Prepare neighbours structure based on cap layout
cd ('/Users/mattiaipem/Documents/IPEM/Projects/ Lou _ Entrainment/Event_Related_Instability');

cfg = [];
cfg.elec = ft_read_sens('standard_waveguard64.elc'); % read electrodes from document (3D coordinates)
cfg.channel = 'all';
cfg.method = 'triangulation';
cfg.compress = 'yes';
cfg.feedback = 'yes';
neighbours = ft_prepare_neighbours(cfg);


% Initialize structure
interpData = struct([]);
% Loop over conditions
for condi = 1:length(postIcaData)
    
% Interpolation (it's important when it comes to compare across subjects)
cfg = [];
cfg.neighbours = neighbours; %defined in previous step
cfg.method = 'average'; 
cfg.layout = layout; 
cfg.badchannel = badChan{condi};
cfg.senstype = 'EEG';

   % Hack capital letters mismatch (match between template and data structure)
        % NB: if you use the hack, revert the upper case afterwards.
        % postIcaData(condi).label = upper(postIcaData(condi).label);

        % Interpolate and assign to temporary variable
        interpTemp = ft_channelrepair(cfg, postIcaData(condi));
        % Re-sort according to fixed template (CPZ end up always appended at the end of the list, followed by reconstructed channels)
        [~,~,idx_resort] = intersect( interpTemp.elec.label , upper(interpTemp.label) ,'stable' );
        %...first the labels...
        interpTemp.label = interpTemp.label(idx_resort);
        interpTemp.trial{1} = interpTemp.trial{1}(idx_resort,:); 

        % Assign to temporary variable
        % ALREADY DONE ABOVE %interpTemp = ft_channelrepair(cfg, postIcaData(condi));

% Concatenate structures
interpData = [interpData interpTemp];

end
% Final check for channel labels (if not all the same, you get an error)
disp([interpData(1).label, interpData(2).label, interpData(3).label]); 
% Prevent that fields have the same length but <64 channels
if length(interpData(1).label) < 64
    warning('Not enough channels! Fix channels before proceeding.')  
    pause()
end


clear interpTemp
clear cleanData
clear postIcaData


% Checkpoint: clean data
cd (['/Volumes/G-DRIVE USB-C/eegData/ERI/ERI_' num2str(sub) '_eeg']);
% Import data with ICA matrix
save(['ERI_' num2str(sub) '_Clean_metronomes']);

close all
clear

% Data are clean!


