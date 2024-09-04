%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ERI: Behavioral pipeline, 'Master' script for eeg IF analyses.
% Alignament strategy is direclty taken from 'ft_trialfun_ERIstudy.m'
% 
% 
% NB: The perturbation onset does NOT coincide with the midi onset.
% In Phase, the gap is systematic, beacause Bart's algorhitm waits for the
% beat before perturbating; In Tempo, they are variable. My approach
% is valid across perturbation types as implemented, accounting for this issue.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
clc

% Set directory - NOTE: always use folder integrated with eeg data (must be up-to-date)
cd('/Users/mattiaipem/Desktop/ERFA_finalAnalysis')

% Settings
srate    = 1000;   %Hz
interval = 600; % ms
target_frex = 1/(interval/1000); % Hz
bpm = target_frex * 60;
timeall =  465001; % total duration in ms
timevec =  [1:timeall]';
pertwin =  5*srate;  % post-perturbation window edge (absolute)
prewin  = .5*srate;  % pre-perturbation window edge  (absolute)
prepert =  4*srate;  % to shift window backwards to perturbation-free period (absolute)
% Design
nsubs  = 20;
nperts = 20; % if actual number excedes, just keep first 20 perturbations
nconds = 11; 
eventlabels = {'Taps','Beats','Phase +','Phase -','Tempo +','Tempo -','End +','End -'};
nevents = length(eventlabels);
condlabels = {'Listen Metro','Listen Instr1','Listem Instr2', 'Listen Lyrics1','Listen Lyrics2', ...
    'Tap Metro Phase','Tap Metro Tempo','Tap Instr Phase','Tap Instr Tempo','Tap Lyrics Phase','Tap Lyrics Tempo'}; 

% Pre-allocate phase timeseries (create a 'grid' of 0s, to fill with taps
% positions; size in known beforehand and equal for every participant)
[beatPhase,tapPhase] = deal( zeros(timeall, 1) ); 
[beatIfN,beatIfP,tapIfN,tapIfP] = deal( nan(prewin+pertwin,nperts,2) ); % vector for perturbation ground-truths (inst freq of stimuli); 3rd dimension for condition (phase/tempo)
% pre-allocate cell for raw data (will be used for eventual fixes later on in the script; e.g., for Tempo/End flippening)
traw = cell(nsubs,1);

for condi = [6 7]   % Tapping conditions only

for subi = 1:nsubs 
       
% Import data
traw{subi} = csvread(['ERI_' num2str(subi) '_Condition' num2str(condi) '.csv']); 

% Check for false starts
false_thresh = 10000; % safe threshold to separate false start from series of valid perturbations 
% account for more false starts
if sum( diff(traw{subi}(:,1)) > false_thresh ) 
    
    % find index(es)
    false_idx = find( diff(traw{subi}(:,1)) > false_thresh );
    % drop everything before last occurrence
    traw{subi}( 1:false_idx(end) , :) = [];

end

% Find 'time 0'
tzero = traw{subi}(traw{subi}(:,2) == 2,1); % temporary assignment, all beats
tzero = tzero(1);               % keep the first beat only
traw{subi}(:,1) = traw{subi}(:,1) - tzero;  % subtract 'time 0' from raw data

% Assign event onsets based on .csv second column (already aligned) 
% TODO: loop for structure, or cell. At the moment, some will be empty
% depending on condition
taps  = traw{subi}(traw{subi}(:,2) == 1,1);
beats = traw{subi}(traw{subi}(:,2) == 2,1); 
pertP = traw{subi}(traw{subi}(:,2) == 3 | traw{subi}(:,2) == 5, 1); % 3 for phase, 5 for tempo
pertN = traw{subi}(traw{subi}(:,2) == 4 | traw{subi}(:,2) == 6, 1); % 4 for phase, 6 for tempo
% 'stop' is logged for tempo perturbations only
pertP_stop = traw{subi}(traw{subi}(:,2) == 7,1);
pertN_stop = traw{subi}(traw{subi}(:,2) == 8,1);
% Get all perturbations together (for further alignment in eeg data)
pertAll = sort([pertN;pertP]);
% Assess distance between perturbation logs and beat logs.
% Find closest neighbour among beats
[idxP,~] = dsearchn(beats,pertP); 
[idxN,~] = dsearchn(beats,pertN);
[idxAll,~] = dsearchn(beats,pertAll);
% Compute distance between perturbation and closest beat
distP(condi-5,subi,:) = beats(idxP) - pertP; %  Note: idx for condi is shifted from [6 7] to [1 2], to adapt to eeg section
distN(condi-5,subi,:) = beats(idxN) - pertN;
distAll(condi-5,subi,:) = beats(idxAll) - pertAll;

% % Inspect in console
% disp([condlabels{condi} 'Logs distances, per + perturbation type'])
% disp([pertP squeeze(distP(condi-5,subi,:))])
% disp([condlabels{condi} 'Logs distances, per - perturbation type'])
% disp([pertN squeeze(distN(condi-5,subi,:))])

% % Correct gap
pertP = pertP - squeeze(abs(distP(condi-5,subi,:))); 
pertN = pertN - squeeze(abs(distN(condi-5,subi,:))); 
fixedAll= sort([pertN;pertP]); % aligned version of pertAll (just use for further inspection)


% Compute derivatives (match size with onsets vector, for further indexing)
iti = [0 ; diff(taps)]; % inter-tap interval - concatenate 0 to match taps size, for debouncing
ibi = [diff(beats)]; % inter-beat interval

% Debounce tapping onsets (false positives a.k.a. double-taps based on ITIs)
% Exclude false positives from onsets timeseries
false_positive = 350;
taps = taps(iti > false_positive | iti == 0); % second condition preserves the 1st tap
% Exclude false positives from ioi timeseries
iti = iti(iti > false_positive | iti ~= 0);  % here we remove the 1st iti (would bias eventual outcome measures)


% % Quality check 
% % Beat standard deviation in ms
% clc
% disp(['Beat mean = ' num2str(mean(ibi)) ' ; standard deviation = ' num2str(std(ibi)) ' ms']);
% 
% % Inspect
% figure;
% subplot(211)
% plot(ibi)
% hold on
% plot(diff(ibi) + mean(ibi))
% title('Inter-beats intervals')
% subplot(212)
% plot(iti)
% title('Inter-taps intervals')
% sgtitle(['Subject: ' num2str(subi) ' - Contition: ' condlabels{condi}])


% Interpolate onsets (do not bother about extremes)
% Exclude eventual negative taps (i.e., taps before time = 0): they are invalid indexes 
taps = taps(taps > 0);
% Put taps in position on the 0s grid, for interpolation
tapPhase(taps(1:end),1) = 1; % no need to skip first tap, as we removed negative values 
% NaN padding (we use non-numerical values, not to bias the results) 
tapPhase(1:taps(1)-1) = nan;      % sacrifice samples up to first tap
tapPhase(taps(end):end) = nan;    % sacrifice samples from last tap until the end

% Compute time series
for i = 1:length(taps)-1

    % Linearly interpolate segments with ramp waves from 0 to 2pi
    tapPhase(taps(i):taps(i+1)) = linspace( 0 , 2*pi , length(taps(i):taps(i+1)));

end


% Repeat for beats
beatPhase(beats(2:end),1) = 1;  % skip 1st beat ('t0', invalid index)
beatPhase(1:beats(1)-1)   = nan;     
beatPhase(beats(end):end) = nan;    

% Compute time series
for i = 2:length(beats)-1

    % Linearly interpolate segments with ramp waves from 0 to 2pi
    beatPhase(beats(i):beats(i+1)) = linspace( 0 , 2*pi , length(beats(i):beats(i+1)));

end


% Compute instantaneous frequencies
tapIf  = [0 ; srate*diff(unwrap(tapPhase))/(2*pi)];
beatIf = [0 ; srate*diff(unwrap(beatPhase))/(2*pi)];


% Assign segments within windows of interest
for perti = 1:nperts-1
    
    beatIfN(:,perti,condi-5) = beatIf( pertN(perti)-prewin:pertN(perti)+pertwin-1 );
    beatIfP(:,perti,condi-5) = beatIf( pertP(perti)-prewin:pertP(perti)+pertwin-1 );
    tapIfN(:,perti,condi-5) = tapIf( pertN(perti)-prewin:pertN(perti)+pertwin-1 );
    tapIfP(:,perti,condi-5) = tapIf( pertP(perti)-prewin:pertP(perti)+pertwin-1 );
    % Perturbation-free events, away from perturbation
    beatIfB(:,perti,condi-5) = beatIf( pertN(perti)-prewin-prepert:pertN(perti)-1 );
    tapIfB(:,perti,condi-5)  =  tapIf( pertN(perti)-prewin-prepert:pertN(perti)-1 );

end

% Perform mean normalization (relative change, in percentage units) 
beatIfN(:,:,condi-5)  = 100 * ( beatIfN(:,:,condi-5)  - mean(beatIfN(1:prewin,:,condi-5),1) ) / target_frex;
beatIfB(:,:,condi-5)  = 100 * ( beatIfB(:,:,condi-5)  - mean(beatIfB(1:prewin,:,condi-5),1) ) / target_frex;
beatIfP(:,:,condi-5)  = 100 * ( beatIfP(:,:,condi-5)  - mean(beatIfP(1:prewin,:,condi-5),1) ) / target_frex;
tapIfN(:,:,condi-5)   = 100 * ( tapIfN(:,:,condi-5)  - mean(tapIfN(1:prewin,:,condi-5),1) ) / target_frex;
tapIfB(:,:,condi-5)   = 100 * ( tapIfB(:,:,condi-5)  - mean(tapIfB(1:prewin,:,condi-5),1) ) / target_frex;
tapIfP(:,:,condi-5)   = 100 * ( tapIfP(:,:,condi-5)  - mean(tapIfP(1:prewin,:,condi-5),1) ) / target_frex;

% Compute and assign average per subject
avgBehavERI_mins{condi-5}(:,subi) = mean(tapIfN(:,:,condi-5),2,'omitnan');
avgBehavERI_null{condi-5}(:,subi) = mean(tapIfB(:,:,condi-5),2,'omitnan');
avgBehavERI_plus{condi-5}(:,subi) = mean(tapIfP(:,:,condi-5),2,'omitnan');

% % Visualize
% figure(condi-5); 
% subplot(4,5,subi) ,
% hold on
% plot(mean(beatIfN,2,'omitnan') , 'r:' , 'LineWidth' , 1.6)
% plot(mean(beatIfP,2,'omitnan') , 'g:' , 'LineWidth' , 1.6)
% plot(mean(beatIfB,2,'omitnan') , 'k:' , 'LineWidth' , 1.6)
% plot(avgIfN(:,subi) , 'r', 'LineWidth' , 1.6)
% plot(avgIfP(:,subi) , 'g', 'LineWidth' , 1.6)
% % plot(mean(tapIfB,2,'omitnan') , 'k', 'LineWidth' , 1.6)
% xline(prewin);
% if condi == 7 || condi == 8
% xline(prewin+3000);
% end
% xlim([0 pertwin])
% xticks([0:500:pertwin])
% xticklabels([0:500:pertwin] - prewin)
% ylim([period-.4*period period+.4*period])
% xlabel('Time (ms)')
% ylabel('Frequency (Hz)')
% % subplot(212) , hold on
% % plot(std(beatIfN,0,2,'omitnan') , 'r:' , 'LineWidth' , 1.6)
% % plot(std(beatIfP,0,2,'omitnan') , 'g:' , 'LineWidth' , 1.6)
% % plot(std(tapIfN,0,2,'omitnan') , 'r', 'LineWidth' , 1.6)
% % plot(std(tapIfP,0,2,'omitnan') , 'g', 'LineWidth' , 1.6)
% title(['Sub #' num2str(subi)])
% sgtitle(condlabels{condi})


% 
% % Save figures for presentation
% % If u do, re-set cd to subjects' folders outside of
% % this loop
% cd('/Users/mattiaipem/Desktop')
% 
% h = figure(1); clf 
% hold on
% plot(mean(beatIfN,2,'omitnan') , 'r:' , 'LineWidth' , 2)
% plot(mean(beatIfP,2,'omitnan') , 'g:' , 'LineWidth' , 2)
% plot(mean(beatIfB,2,'omitnan') , 'k:' , 'LineWidth' , 2)
% plot(avgIfN(:,subi) , 'r', 'LineWidth' , 2)
% plot(avgIfP(:,subi) , 'g', 'LineWidth' , 2)
% % plot(mean(tapIfB,2,'omitnan') , 'k', 'LineWidth' , 1.6)
% xline(prewin);
% if condi == 7 || condi == 8
% xline(prewin+3000);
% end
% xlim([0 pertwin])
% xticks([0:500:pertwin])
% xticklabels([0:500:pertwin] - prewin)
% ylim([period-.4*period period+.4*period])
% xlabel('Time (ms)')
% ylabel('Frequency (Hz)')
% title(['Sub #' num2str(subi)])
% axis square
% saveas(h,[condlabels{condi} '_sub#' num2str(subi) '.jpg'])





% % Inspection #2
% figure;
% plot(beatIf)
% hold on
% for i = 1:length(pertAll)
%    xline(pertAll(i)); 
%    xline(fixedAll(i), 'r');
% end
% ylim([target_frex-.2*target_frex target_frex+.2*target_frex])
% title(['Misalignment check: ' condlabels{condi} ' - Sub #' num2str(subi)])
% 
% 
% 

end

end


% Plot Grand-Average with MSE
fontsize = 25;
linethicker = 3;
onsetthicker = 1.5;
% vectors for tempo + and -
[tminsvec,tplusvec] = deal( zeros(length(avgBehavERI_mins{1}),1) );
tminsvec(prewin:prewin+3000) = -10;
tplusvec(prewin:prewin+3000) = 10;

for condi = [6 7] % 6 = phase shift; 7 = tempo change

h_gravg = figure(100+condi) ; clf,  hold on
errorbar(mean(avgBehavERI_mins{condi-5},2) , std(avgBehavERI_mins{condi-5},0,2,'omitnan')/sqrt(nsubs-1) , 'r' , 'LineWidth' , .5)
errorbar(mean(avgBehavERI_plus{condi-5},2) , std(avgBehavERI_plus{condi-5},0,2,'omitnan')/sqrt(nsubs-1) , 'g' , 'LineWidth' , .5)
plot(mean(avgBehavERI_mins{condi-5},2) , 'k' , 'LineWidth' , linethicker)
plot(mean(avgBehavERI_null{condi-5},2) , 'k' , 'LineWidth' , linethicker)
plot(mean(avgBehavERI_plus{condi-5},2) , 'k' , 'LineWidth' , linethicker)
% plot(mean(tapIfB,2,'omitnan') , 'k', 'LineWidth' , 1.6)
xline(prewin);
yline(0);
if condi == 7 || condi == 8
     plot(tminsvec,'r--', 'LineWidth', linethicker)
     plot(tplusvec,'g--', 'LineWidth', linethicker)
     ylim([-20 20])
else
    plot( mean( (beatIfN(:,:,condi-5)-mean(max(beatIfN(:,:,condi-5)),'omitnan')) ,2,'omitnan' ) , 'r--' , 'LineWidth' , linethicker) % correct slight offset in groundtruth, must be a very minor mistake
    plot( mean( beatIfP(:,:,condi-5),2,'omitnan' ) , 'g--' , 'LineWidth' , linethicker)
    ylim([-40 40])
end
xlim([0 pertwin])
xticks([500:1000:pertwin]) 
xticklabels([500:1000:pertwin] - prewin)
xlabel('Time (ms)')
ylabel('Frequency (% change)')
axis square
ax = gca;
ax.LineWidth = 3;
set(ax,'fontsize', fontsize)

%saveas(h_gravg,[condlabels{condi} '_Grand-Average.jpg'])

end







%clearvars -except beatIfN beatIfB beatIfP tapIfN tapIfB tapIfP avgBehavERI_mins avgBehavERI_null avgBehavERI_plus distN distP distAll prepert traw eventlabels






%% GED (process entrained component in Fieldtrip structure)

% From now on, we only process Tap conditions and use Listen to compute
% auditory entrained component. Produce 2 outputs (gedMotorData,ged
% AudioData) at the end of the block. Trials ERI will be analysed
% separately for the 2 components.


close all

% Re-initialize some settings
nsubs  = 20;
nchans = 64;
srate  = 1000;

% Define GED configuration
cfgGed = [];
cfgGed.targetfrex = 1.667; % Beat frequency (100 BPM)
cfgGed.fwhm       = .3; % full width at hald maximum (allows for frequency fluctuations)
% Set the following to zero if R is not computed from neighbours
cfgGed.neig      = 0;  % distance to frequency neighbors in Hz
cfgGed.fwhm_neig = 0;  % fwhm in Hz for neighbors
cfgGed.shr    = .01;      % shrinkage proportion
cfgGed.srate  = srate;  %sampling rate
cfgGed.winlen = 60;    % time window length, in seconds (steady period)  
cfgGed.covwin = [0 .6*cfgGed.srate];     % time window for covarince matrix 
cfgGed.chans2keep = 'all';  % by nature of the experiment, there should be no preconception about ROIs (data-drive separation of Audio-Motor components)
cfgGed.comp2pick  = 1;       % select component based on the ranking of eigenvalues 


% Initializer GED outputs, to store per participant
[gedAudioMap,gedMotorMap] = deal(nan(nsubs,nchans)); % use NaNs, and 'omitnan' for grand-average (account for lost participants)
[gedAudioSNR,gedMotorSNR] = deal(nan(cfgGed.winlen*srate,nsubs));
[gedAudioEvals,gedMotorEvals] = deal(nan(nchans,nsubs));

% Pick subject / loop over subjects
for subi = [1:5 7:20]


% Import data 
load(['ERI_' num2str(subi) '_Clean_metronomes']);


% Initialize FT structure to store GED outputs
gedMotorData = interpData; % 'FieldTrick': copy Data and replace channels; ft format needed to redefine trial later 
gedAudioData = interpData; 


% Define "Listen" condition (for Metronomes, this is 1; for Music, depends on the mix)
lstncond = 1; 

% Loop over tapping conditions
for condi = 2:3 %tapconds % we only process Tap conditions
    
    % Assign GED input: multivariate timeseries (macro-trial)
    input_ged = interpData(condi).trial{1};    
  
    % Compute GED output #1: Sensorimotor component
    [gedMotorCompt,gedMotorComptX,gedMotorSNR(:,subi),gedMotorMap(subi,:),gedMotorEvals(:,subi)] = ft_myGed(cfgGed , input_ged , 0);
    
    % Concatenate stable period of Listen, to use as the first cfg.timewin seconds 
    input_ged = [interpData(lstncond).trial{1} interpData(condi).trial{1}]; % it is crucial that channels are the same across datasets 
    % TODO: for Music conditions, deal wiht randomization of mixes 1&2 (ideally, this should be done early in the pre-processing, or even before the very beginning of it); would it be the case, delete this comment when the time comes    
    % Compute GED output #2: Auditory component (NB: input includes period of Listening in the beginning)
    [gedAudioCompt,gedAudioComptX,gedAudioSNR(:,subi),gedAudioMap(subi,:),gedAudioEvals(:,subi)] = ft_myGed(cfgGed , input_ged , 0);
    % Drop stable period of Listen in the beginning
    win2drop = [1:length(interpData(lstncond).trial{1})];
    gedAudioCompt(win2drop) = [];
    %NB: spectra are still not matching in size! You can't just drop
    %them... the reason being that they are computed inside ft_myGed()
    %starting from timeseries of different length (due to concatenation of Listening Condition)
    % I considered swapping the 1st minute (Tapping for Listening) instead
    % of concatenating and dropping... but it messes things up. There must
    % have been a reason why I did not do it before
     
    % Assign timeseries, replacing interpData (see 'Fieldtrick' above)
    gedAudioData(condi).trial{1} = gedAudioCompt;    
    gedMotorData(condi).trial{1} = gedMotorCompt;  
    % Remove channel labels, replace with single component lable 'Audio', 'Motor'
    gedAudioData(condi).label = {'Audio'}; 
    gedMotorData(condi).label = {'Motor'};


end

% NOTE: from now on, structures will only contain "Tap" conditions
gedAudioData(1) = [];
gedMotorData(1) = [];

%% Bad trial rejections on component timeseries

% 'FLIPPENING' ISSUE: we check and correct for that here

% Check based on double logic consition
logslogic = [nan diff([cfgTrial(3).event.sample])] < 3050; % pick Tempo condition; leave a 50ms margin
labslogic = contains({cfgTrial(3).event.value},'End');     
fixlogic  = logslogic' - labslogic';                    % where incongruency happens
if sum(abs(fixlogic)) ~= 0
    warning(' "Flippening" issue happened! Wrong End-Tempo pairs have been reverted.')
       
    % replace events with Teensy logs (these are always correct); 
    % defined above in eventlabels cell    
    tmpend = traw{subi}(traw{subi}(:,2)>4 ,1:2) % contains timestamps and labels of Tempo and End triggers
    tmpend(2:end,3) = diff(tmpend(:,1)); % add inter-triggers interval
    % figure
    % plot(diff(tmpend(:,2)))
    % plot(diff(tmpend(:,1)))

    % find where flippening happens
    flipidx = sum(fixlogic==0) + 1;
    
% Inspect    
figure(300+subi);
subplot(211)
plot( diff([cfgTrial(3).event.sample]) )
xline(flipidx(1),'r:')
{cfgTrial(3).event.type}
[cfgTrial(3).event.value]
xticklabels([{cfgTrial(2).event.value}])
xticklabels({})
xticks(1:length(cfgTrial(3).event))
xticklabels({cfgTrial(3).event.value})
xtickangle(90)
title('Before reverting flipped pairs')

    % from flippening onwards,
    for flipi = flipidx:length(fixlogic)
        cfgTrial(3).event(flipi).value = eventlabels{tmpend(flipi,2)}; % replace flipped labels, accounting for gap in vector sizes (event structure misses 1 element)
    end

subplot(212)
plot( diff([cfgTrial(3).event.sample]) )
xline(flipidx(1),'r:')
{cfgTrial(3).event.type}
[cfgTrial(3).event.value]
xticklabels([{cfgTrial(3).event.value}])
xticklabels({})
xticks(1:length(cfgTrial(3).event))
xticklabels({cfgTrial(3).event.value})
xtickangle(90)
title('After reverting flipped pairs')

  
else
    disp('No issue with Tempo/End flippening')
    
end
% NOTE: we are not inserting the missing 'End', as we are ignoring such
% events and only basing our segmentation on 'Tempo'


% NB store bad trials, but REMOVE ONLY AT THE VERY VERY END OF PROCESSING
% (avoid segmentation for inst. frex processing)

% Drop 'Listen' conditions from trial configuration
cfgTrial(1:length(lstncond)) = [];

% Pre-allocate ERP structures
% [gedAudioERP, gedMotorERP] = deal(struct([]));

% Loop over conditions (in every block). Note 'Listening' is out of the
% count now
for condi = 1:2
        
    % Re-define trials for ERP
    % Update cfgTrial structure
    cfgTrial(condi).tzero = cfgTrial(condi).event(1).sample; % assign 'time 0' ('Start' sample)
    cfgTrial(condi).event(1) = []; % drop 'Start'
    % Drop 'End' events, in order to make trial removal more
    % straightforward (would you encounter issues with trial removal, there
    % might be a mistake here to double check)
    cfgTrial(condi).event(contains({cfgTrial(condi).event.value},'End')) = [];
    cfgTrial(condi).erwinsize = 5 * srate; % event-related windowsize
    % Assign samples of start and end of trials; 3rd entry is 'offset' (=0)
    cfgTrial(condi).trl = [[cfgTrial(condi).event.sample]'  - cfgTrial(condi).tzero , [cfgTrial(condi).event.sample]' - cfgTrial(condi).tzero  + cfgTrial(condi).erwinsize , zeros(length(cfgTrial(condi).event),1)]; 
    cfgTrial(condi).trialdef.eventtype  = 'Stimulus';
    cfgTrial(condi).trialdef.eventvalue = {'Phase -' , 'Phase +' , 'Tempo -' , 'Tempo +'}; % Exclude 'End' events, not to mess up the trial removal (otherwise you will need piar-wise removal in Tempo condition)
    % Correct misalignment perturbation-closest beat; based on distances
    % computed in first behavioral block
    cfgTrial(condi).trl(:,1:2) = cfgTrial(condi).trl(:,1:2) - squeeze(distAll(condi,subi,:));
    %cfgTrial.trialfun = 'ft_trialfun_general';
    %cfgTrial.trl = cfgTrial.trl - cfgTrial.tzero % I think I'd better not
    %re-align to time zero, but rather keet the original timestamps for
    %alignment with behavior
    % Segmenting function with configuration and continuous filtered data as input argument
    gedAudioERP(condi) = ft_redefinetrial(cfgTrial(condi) , gedAudioData(condi)); 
    gedMotorERP(condi) = ft_redefinetrial(cfgTrial(condi) , gedMotorData(condi));
    
    % Define perturbation-free trials, for baseline (NOTE: we are overwriting the variable; bad practice, but it works so far)
    cfgTrial(condi).trl(:,1:2) = cfgTrial(condi).trl(:,1:2) - prepert; % shift onsets
    % Segmenting function with configuration and continuous filtered data as input argument
    gedAudioERP_Bl(condi) = ft_redefinetrial(cfgTrial(condi) , gedAudioData(condi));
    gedMotorERP_Bl(condi) = ft_redefinetrial(cfgTrial(condi) , gedMotorData(condi)); 
              
%     % Refresh cfg
%     cfg = [];
%     cfg.method = 'summary'; %(?)
%     cfg.alim = 5e-5;     %set limits for the graph
%     gedAudioERP(condi) = ft_rejectvisual(cfg,gedAudioERP(condi));
%     pause;
%     gedMotorERP(condi) = ft_rejectvisual(cfg,gedMotorERP(condi));
%     pause;


end

%% Insert bad trials

% Bad trials are not stored automatically; store manually but DO NOT REMOVE
% YET.
% % NB: Put phase perturbations (left entries first!)
                  

badAudioTrial  = {[]         []};

badMotorTrial  = {[]         []};


% badAudioTrial  = {[]         []};
%                                  
% badMotorTrial  = {[]         []};


% badAudioTrial  = {[]         []};
%                                  
% badMotorTrial  = {[]         []};


%% Instantaneous frequency computation

% Compute on macrotrial; remove afterwards the bad trials previously marked

% 'Fieldtrick': mimick FT structure
ifMotorData = gedMotorData; 
ifAudioData = gedAudioData; 

% Initialize structure
%instfrexData = struct([]);

% Configuration for Instantaneous FGrequency
cfgIf = [];
cfgIf.targetfrex = cfgGed.targetfrex; % same as GED filter
cfgIf.fwhm = cfgGed.fwhm;
% Set steps for median smoothing
cfgIf.medsamples = 400; % 400ms recommended by Cohen(2014); also order of 10

% Compute instantaneous frequency
for condi = 1:2     

% First on component...    
% Re-filter the component (choose your filter)
tempIf = filterFGx( gedMotorData(condi).trial{1}, srate, cfgIf.targetfrex, cfgIf.fwhm, 0); %Gaussian filter
% Extract angle time series
tempIf = angle(hilbert(tempIf));
% Compute instantaneous frequency
tempIf = [0 srate*diff(unwrap(tempIf))/(2*pi)];
% Median filter (remove the abrupt peaks, which are normal)
tempIf = movmedian(tempIf,cfgIf.medsamples);

% replace long timeseries        
ifMotorData(condi).trial{1} = tempIf;


% ... then repeat for Audio component
tempIf = filterFGx( gedAudioData(condi).trial{1}, srate, cfgIf.targetfrex, cfgIf.fwhm, 0); %Gaussian filter
tempIf = angle(hilbert(tempIf));
tempIf = [0 srate*diff(unwrap(tempIf))/(2*pi)];
tempIf = movmedian(tempIf,cfgIf.medsamples);
ifAudioData(condi).trial{1} = tempIf;


end

clear tempIf

%% Bad trials removal

%NB: DO NOT RUN THIS BLOCK TWICE, BECAUSE YOU'D KEEP DELETING TRIALS
%EVERYTIME. SOLVE THIS

% NOTE2: when removing trials, remove from distaAll accordingly for
% aligning perturbation-event


% Based on previous identification in ERP (variance in voltage as criterion for removal)

% Define event labels, for different conditions
pluslabels = {'Phase +' , 'Tempo +'}; 
minslabels = {'Phase -' , 'Tempo -'};

for condi = 1:2

    % NB: the following trial definition needs to be the exact SAME AS ERP. 
    % check that nothing changed in that previous block

    % Re-define trials ERI
    % Update cfgTrial structure
    cfgTrial(condi).tzero = 0; %cfgTrial.event(1).sample; % assign 'time 0' ('Start' sample)
    cfgTrial(condi).event(1) = []; % drop 'Start'
    cfgTrial(condi).prestim = .5*srate;
    cfgTrial(condi).erwinsize = 5 * srate; % event-related windowsize
    % Assign samples of start and end of trials; 3rd entry is 'offset' (=0)
    cfgTrial(condi).trl = [[cfgTrial(condi).event.sample]'  - cfgTrial(condi).prestim , [cfgTrial(condi).event.sample]'  + cfgTrial(condi).erwinsize , zeros(length(cfgTrial(condi).event),1)]; 
    % Correct misalignment perturbation-closest beat; based on distances
    % computed in first behavioral block
    cfgTrial(condi).trl(:,1:2) = squeeze(distAll(condi,subi,1:length(cfgTrial(condi).trl))) + cfgTrial(condi).trl(:,1:2);
    %cfgTrial.trialfun = 'ft_trialfun_general';
    %cfgTrial.trl = cfgTrial.trl - cfgTrial.tzero
    % Segmenting function with Trial configuration and continuous filtered data as input argument
    gedMotorERI(condi) = ft_redefinetrial(cfgTrial(condi) , ifMotorData(condi)); 
    gedAudioERI(condi) = ft_redefinetrial(cfgTrial(condi) , ifAudioData(condi));     
    % Define perturbation-free trials, for baseline
    cfgTrial(condi).trl(:,1:2) = cfgTrial(condi).trl(:,1:2) - prepert; % shift onsets
    % Segmenting function with Baseline configuration and continuous filtered data as input argument
    gedMotorERI_Bl(condi) = ft_redefinetrial(cfgTrial(condi) , ifMotorData(condi));
    gedAudioERI_Bl(condi) = ft_redefinetrial(cfgTrial(condi) , ifAudioData(condi)); 
    
    
    % Remove bad trials
    cfgERI = [];
    cfgERI.vector = [1:length(gedMotorERI(condi).trial)];
    cfgERI.trials = cfgERI.vector(~ismember(cfgERI.vector,badMotorTrial{condi}));  
    gedMotorERI(condi) = ft_preprocessing(cfgERI , gedMotorERI(condi));
    
    % Repeat for audio component
    cfgERI = [];
    cfgERI.vector = [1:length(gedAudioERI(condi).trial)];
    cfgERI.trials = cfgERI.vector(~ismember(cfgERI.vector,badAudioTrial{condi}));  
    gedAudioERI(condi) = ft_preprocessing(cfgERI , gedAudioERI(condi));
    

    % WATCH OUT FOR MISALIGNMENT!
    % Visualise
    % cfg = [];
    % cfg.ylim = [cfgGed.targetfrex - .5 , cfgGed.targetfrex + .5];
    % ft_databrowser(cfg,gedMotorERI(1))


% Get out of FT format, select ERIs per type 
% Find indexes
idxPlus = find(ismember({cfgTrial(condi).event.value}, pluslabels{condi}));
idxMins = find(ismember({cfgTrial(condi).event.value}, minslabels{condi}));


% Computation and assignment of ERPs and ERIs
for triali = 1:floor(ntrials/2) - 1 % - length(badChan{condi})
    % Assign 
    % ERPs
    gedMotorERP_plus{condi}(triali,:) = gedMotorERP(condi).trial{ idxPlus(triali) };
    gedMotorERP_mins{condi}(triali,:) = gedMotorERP(condi).trial{ idxMins(triali) };
    gedAudioERP_plus{condi}(triali,:) = gedAudioERP(condi).trial{ idxPlus(triali) };
    gedAudioERP_mins{condi}(triali,:) = gedAudioERP(condi).trial{ idxMins(triali) };
    % ERIs
    gedMotorERI_plus{condi}(triali,:) = gedMotorERI(condi).trial{ idxPlus(triali) };
    gedMotorERI_mins{condi}(triali,:) = gedMotorERI(condi).trial{ idxMins(triali) };
    gedAudioERI_plus{condi}(triali,:) = gedAudioERI(condi).trial{ idxPlus(triali) };
    gedAudioERI_mins{condi}(triali,:) = gedAudioERI(condi).trial{ idxMins(triali) };
end
% now for perturbation-free 'null' trials
for triali = 1:(floor(ntrials/2) - 1) * 2 % the amount of steady segments is double of individual directions
    % ERPs
    gedMotorERP_null{condi}(triali,:) = gedMotorERP_Bl(condi).trial{triali};
    gedAudioERP_null{condi}(triali,:) = gedAudioERP_Bl(condi).trial{triali};
    % ERIs
    gedMotorERI_null{condi}(triali,:) = gedMotorERI_Bl(condi).trial{triali};
    gedAudioERI_null{condi}(triali,:) = gedAudioERI_Bl(condi).trial{triali};
end
    
    
% Baseline subtraction/normalization 
% ERPs
gedMotorERP_plus{condi} = gedMotorERP_plus{condi} - mean ( gedMotorERP_plus{condi}(:,1:cfgTrial(condi).prestim) , 2 );
gedMotorERP_mins{condi} = gedMotorERP_mins{condi} - mean ( gedMotorERP_mins{condi}(:,1:cfgTrial(condi).prestim) , 2 );
gedMotorERP_null{condi} = gedMotorERP_null{condi} - mean ( gedMotorERP_null{condi}(:,1:cfgTrial(condi).prestim) , 2 );
gedAudioERP_plus{condi} = gedAudioERP_plus{condi} - mean ( gedAudioERP_plus{condi}(:,1:cfgTrial(condi).prestim) , 2 );
gedAudioERP_mins{condi} = gedAudioERP_mins{condi} - mean ( gedAudioERP_mins{condi}(:,1:cfgTrial(condi).prestim) , 2 );
gedAudioERP_null{condi} = gedAudioERP_null{condi} - mean ( gedAudioERP_null{condi}(:,1:cfgTrial(condi).prestim) , 2 );

% ERIs (same normalization as behavioral data; expressed in percent change units)
gedMotorERI_plus{condi} = 100 * (gedMotorERI_plus{condi} - mean ( gedMotorERI_plus{condi}(:,1:cfgTrial(condi).prestim) , 2 ) ) / cfgGed.targetfrex;
gedMotorERI_mins{condi} = 100 * (gedMotorERI_mins{condi} - mean ( gedMotorERI_mins{condi}(:,1:cfgTrial(condi).prestim) , 2 ) ) / cfgGed.targetfrex;
gedMotorERI_null{condi} = 100 * (gedMotorERI_null{condi} - mean ( gedMotorERI_null{condi}(:,1:cfgTrial(condi).prestim) , 2 ) ) / cfgGed.targetfrex;
gedAudioERI_plus{condi} = 100 * (gedAudioERI_plus{condi} - mean ( gedAudioERI_plus{condi}(:,1:cfgTrial(condi).prestim) , 2 ) ) / cfgGed.targetfrex;
gedAudioERI_mins{condi} = 100 * (gedAudioERI_mins{condi} - mean ( gedAudioERI_mins{condi}(:,1:cfgTrial(condi).prestim) , 2 ) ) / cfgGed.targetfrex;
gedAudioERI_null{condi} = 100 * (gedAudioERI_null{condi} - mean ( gedAudioERI_null{condi}(:,1:cfgTrial(condi).prestim) , 2 ) ) / cfgGed.targetfrex;


% Compute average within-subjects
% ERPs
avgMotorERP_plus{condi}(subi,:) = mean(gedMotorERP_plus{condi},1);
avgMotorERP_mins{condi}(subi,:) = mean(gedMotorERP_mins{condi},1);
avgMotorERP_null{condi}(subi,:) = mean(gedMotorERP_null{condi},1);
avgAudioERP_plus{condi}(subi,:) = mean(gedAudioERP_plus{condi},1);
avgAudioERP_mins{condi}(subi,:) = mean(gedAudioERP_mins{condi},1);
avgAudioERP_null{condi}(subi,:) = mean(gedAudioERP_null{condi},1);

% ERIs
avgMotorERI_plus{condi}(subi,:) = mean(gedMotorERI_plus{condi},1,'omitnan');
avgMotorERI_mins{condi}(subi,:) = mean(gedMotorERI_mins{condi},1,'omitnan');
avgMotorERI_null{condi}(subi,:) = mean(gedMotorERI_null{condi},1,'omitnan');
avgAudioERI_plus{condi}(subi,:) = mean(gedAudioERI_plus{condi},1,'omitnan');
avgAudioERI_mins{condi}(subi,:) = mean(gedAudioERI_mins{condi},1,'omitnan');
avgAudioERI_null{condi}(subi,:) = mean(gedAudioERI_null{condi},1,'omitnan');





% Visualization 

% Visualize Motor component
eriTime = [-cfgTrial(condi).prestim:cfgTrial(condi).erwinsize+cfgTrial(condi).prestim];
h_motor = figure(40+condi) ; clf , hold on
% Instantaneous frequency response
% Single trials
for triali = 1:size(gedAudioERI_plus{condi},1)

plot(gedAudioERI_plus{condi}(triali,:),'g:','LineWidth',1.5)
plot(gedAudioERI_mins{condi}(triali,:),'r:','LineWidth',1.5)

end
% Average trial
plot(avgMotorERI_plus{condi}(subi,:) ,'g','LineWidth',2.5  )
plot(avgMotorERI_mins{condi}(subi,:) ,'r','LineWidth',2.5  )
plot(avgMotorERI_null{condi}(subi,:) ,'o','LineWidth',2.5  )
yline(0,'k:');
xticks([0:500:length(avgMotorERI_mins{condi})-1])
xticklabels(num2str([-cfgTrial(condi).prestim:500:5500]'))
xline(cfgTrial(condi).prestim);
if condi == 7 || condi == 8
xline(prewin+3000);
end
xlabel('Time(ms)')
ylabel('Frequency (% change)')
% ylim([-.2 .2])
axis square
legend({pluslabels{condi},minslabels{condi}})
title(['Sub #' num2str(subi)])

% Visualize Audio component
eriTime = [-cfgTrial(condi).prestim:cfgTrial(condi).erwinsize+cfgTrial(condi).prestim];
h_audio = figure(50+condi) ; clf , hold on
% Instantaneous frequency response
% Single trials
for triali = 1:size(gedMotorERI_plus{condi},1)

plot(gedMotorERI_plus{condi}(triali,:),'g:','LineWidth',1.5)
plot(gedMotorERI_mins{condi}(triali,:),'r:','LineWidth',1.5)

end% Average trial
plot(avgAudioERI_plus{condi}(subi,:) ,'g','LineWidth',2.5  )
plot(avgAudioERI_mins{condi}(subi,:) ,'r','LineWidth',2.5  )
plot(avgAudioERI_null{condi}(subi,:) ,'o','LineWidth',2.5  )
yline(0,'k:');
xticks([0:500:length(gedAudioERI_mins{condi})-1])
xticklabels(num2str([-cfgTrial(condi).prestim:500:5500]'))
xline(cfgTrial(condi).prestim);
if condi == 7 || condi == 8
xline(prewin+3000);
end
xlabel('Time(ms)')
ylabel('Frequency (% change)')
% ylim([-.2 .2])
axis square
legend({pluslabels{condi},minslabels{condi}})
title(['Sub #' num2str(subi)])


% Save figures
% saveas(h_motor,['motorERI_' condlabels{condi} '_sub#' num2str(sub) '.jpg'])
% saveas(h_audio,['audioERI_' condlabels{condi} '_sub#' num2str(sub) '.jpg'])
% 


end

end

%%
% Plot Grand-Average ERFA with MSE
% Visualize
nsubs = 20;

cond2plot = 2; % 1 = phase shift; 2 = tempo change
fontsize = 25;
linethicker = 3;
onsetthicker = 1.5;
% vectors for tempo + and -
[tminsvec,tplusvec] = deal( zeros(length(avgMotorERI_mins{2}),1) );
tminsvec(cfgTrial(2).prestim:cfgTrial(2).prestim+3000) = -10;
tplusvec(cfgTrial(2).prestim:cfgTrial(2).prestim+3000) = 10;
% Motor component
for condi = cond2plot
    
    figure(100+condi) ; clf,  hold on
    errorbar(mean(avgMotorERI_mins{condi},1) , std(avgMotorERI_mins{condi},0,1,'omitnan')/sqrt(nsubs-1) , 'r' , 'LineWidth' , .5)
    errorbar(mean(avgMotorERI_plus{condi},1) , std(avgMotorERI_plus{condi},0,1,'omitnan')/sqrt(nsubs-1) , 'g' , 'LineWidth' , .5)
    %errorbar(mean(avgMotorERI_null{condi},1) , std(avgMotorERI_null{condi},0,1,'omitnan')/sqrt(nsubs-1) , 'k' , 'LineWidth' , .5)   
    plot(mean(avgMotorERI_mins{condi},1) , 'k' , 'LineWidth' , linethicker)
    plot(mean(avgMotorERI_plus{condi},1) , 'k' , 'LineWidth' , linethicker)
    plot(mean(avgMotorERI_null{condi},1) , 'k' , 'LineWidth' , linethicker)
    ylim([-12 12])
    % plot(mean(tapIfB,2,'omitnan') , 'k', 'LineWidth' , 1.6)
    xlim([0 length(gedAudioERI_mins{condi})-101])
    xticks([500:1000:length(gedAudioERI_mins{condi})-101])
    xticklabels(num2str([-cfgTrial(condi).prestim+500:1000:4500]')) 
    
    xline(cfgTrial(condi).prestim,'LineWidth', onsetthicker);
    if condi >1
         plot(tminsvec,'r--', 'LineWidth', linethicker)
         plot(tplusvec,'g--', 'LineWidth', linethicker)
    else
        plot( mean( (beatIfN(:,:,condi)-mean(max(beatIfN(:,:,condi)),'omitnan')) ,2,'omitnan' ) , 'r--' , 'LineWidth' , linethicker) % correct slight offset in groundtruth, must be a very minor mistake
        plot( mean( beatIfP(:,:,condi),2,'omitnan' ) , 'g--' , 'LineWidth' , linethicker)
    end
    xlabel('Time (ms)')
    ylabel('Frequency (% change)')
    ylim([-12 12])
    axis square
    ax = gca;
    ax.LineWidth = 3
    set(ax,'fontsize', fontsize)
    %legend({pluslabels{condi},minslabels{condi}})
    %     saveas(h_gravg,[condlabels{condi} '_Grand-Average.jpg'])

end

% Auditory component
for condi = cond2plot
    
    figure(200+condi) ; clf,  hold on
    errorbar(mean(avgAudioERI_mins{condi},1) , std(avgAudioERI_mins{condi},0,1,'omitnan')/sqrt(nsubs-1) , 'r' , 'LineWidth' , .5)
    errorbar(mean(avgAudioERI_plus{condi},1) , std(avgAudioERI_plus{condi},0,1,'omitnan')/sqrt(nsubs-1) , 'g' , 'LineWidth' , .5)
%     errorbar(mean(avgAudioERI_null{condi},1) , std(avgAudioERI_null{condi},0,1,'omitnan')/sqrt(nsubs-1) , 'o' , 'LineWidth' , .5)   
    plot(mean(avgAudioERI_mins{condi},1) , 'k' , 'LineWidth' , linethicker)
    plot(mean(avgAudioERI_plus{condi},1) , 'k' , 'LineWidth' , linethicker)
    plot(mean(avgAudioERI_null{condi},1) , 'k' , 'LineWidth' , linethicker)
    ylim([-12 12])
    % plot(mean(tapIfB,2,'omitnan') , 'k', 'LineWidth' , 1.6)
    yline(0,'k:');
    xlim([0 length(gedAudioERI_mins{condi})-101])
    xticks([500:1000:length(gedAudioERI_mins{condi})-101])
    xticklabels(num2str([-cfgTrial(condi).prestim+500:1000:4500]')) 

    xline(cfgTrial(condi).prestim,'LineWidth', onsetthicker);
    if condi >1
         plot(tminsvec,'r--', 'LineWidth', 2)
         plot(tplusvec,'g--', 'LineWidth', 2)
    else
        plot( mean( (beatIfN(:,:,condi)-mean(max(beatIfN(:,:,condi)),'omitnan')) ,2,'omitnan' ) , 'r--' , 'LineWidth' , linethicker) % correct slight offset in groundtruth, must be a very minor mistake
        plot( mean( beatIfP(:,:,condi),2,'omitnan' ) , 'g--' , 'LineWidth' , linethicker)
    end
    xlabel('Time (ms)')
    ylabel('Frequency (% change)')
    ylim([-12 12])
    axis square
    ax = gca;
    ax.LineWidth = 3
    set(ax,'fontsize', fontsize)
%legend({pluslabels{condi},minslabels{condi}})
%     saveas(h_gravg,[condlabels{condi} '_Grand-Average.jpg'])

end

% 
% smth = 100
% % Plot Grand-Average ERPs
% for condi = 1:2
%     
%     figure(300+condi); clf
%     %subplot(211)
%     hold on
%     %plot(smooth(mean(avgAudioERP_mins{condi},1,'omitnan'),smth),'r--')
%     plot(smooth(mean(avgAudioERP_null{condi},1,'omitnan'),smth),'k--')
%     %plot(smooth(mean(avgAudioERP_plus{condi},1,'omitnan'),smth),'g--')
%     %subplot(212)
%     hold on
%     %plot(smooth(mean(avgMotorERP_mins{condi},1,'omitnan'),smth),'r')
%     plot(smooth(mean(avgMotorERP_null{condi},1,'omitnan'),smth),'k')
%     %plot(smooth(mean(avgMotorERP_plus{condi},1,'omitnan'),smth),'g')
%     
% end



%% Plot Grand-Average GED outputs

% Normalize topographical maps
zAudioMap = zscore(gedAudioMap([1:2,4:5,7:20],:),0,2);
zMotorMap = zscore(gedMotorMap([1:2,4:5,7:20],:),0,2);
% prepare frequency vector
pnts = cfgGed.winlen*cfgGed.srate; %length sufficiently high to guarantee resolution at denominator (= pnts)
hz   = linspace(0,cfgGed.srate,pnts); %vector of frequencies; shortcut: all the frex above nyquist are not valid


% Load EEGLAB input for topoplotIndie
fontsize = 18;
load('chanlocs_waveguard64.mat','EEG'); 
figure,clf
%  Audio Topography 
subplot(231)
topoplotIndie(mean(zAudioMap,1) , EEG.chanlocs,'numcontour',0,'electrodes','off');
title('Spatial pattern') 
colorbar
caxis([-.8 .8])
set(gca,'fontsize', fontsize)
%  Audio Spectrum
subplot(232)
plot(hz, smooth( mean(gedAudioSNR,2,'omitnan') , 1),'rs-','markerfacecolor','k','linew',1.5)
set(gca,'xlim',[0 15])     
xlabel('Frequency (Hz)'), ylabel('Power (SNR units)') %or 'Signal-to-noise ratio (%)'
title('SNR spectrum') %or 'SNR' spectrum
ylim([0 170])
axis square
set(gca,'fontsize', fontsize)
%  Audio Eigenspectrum
subplot(233);
plot(squeeze(mean(gedAudioEvals,2,'omitnan')),'rs-','markerfacecolor','k','linew',2)
xlabel('Component'), ylabel('Explained variance (%)')
title('Eigenspectrum')
ylim([0 17])
set(gca,'fontsize', fontsize)
axis square
% Motor Topography
subplot(234)
topoplotIndie(mean(zMotorMap,1) , EEG.chanlocs,'numcontour',0,'electrodes','off');
title('Spatial pattern') 
colorbar
caxis([-.8 .8])
set(gca,'fontsize', fontsize)
% Motor Spectrum
subplot(235)
plot(hz, smooth( mean(gedMotorSNR,2,'omitnan') , 1),'s-','markerfacecolor','k','linew',1.5)
set(gca,'xlim',[0 17])     
xlabel('Frequency (Hz)'), ylabel('Power (SNR units)') %or 'Signal-to-noise ratio (%)'
title('SNR spectrum') %or 'SNR' spectrum
ylim([0 170])
axis square
set(gca,'fontsize', fontsize)
% Motor Eigenspectrum
subplot(236);
plot(squeeze(mean(gedMotorEvals,2,'omitnan')),'s-','markerfacecolor','k','linew',2)
xlabel('Component'), ylabel('Explained variance (%)')
title('Eigenspectrum')
ylim([0 17])
axis square
set(gca,'fontsize', fontsize)



% sgtitle('Grand-average GED outputs')



