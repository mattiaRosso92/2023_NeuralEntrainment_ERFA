% Export for stats in R
clear all
close all
clc

% Import .mat dataset 
load('ERFA_Metronomes_results.mat')

% Calculate integrals for ERPAs ('phase' perturbations only)

% THESE ARE ALREADY INITIALIZED IN MAIN SCRIPT (DOUBLE-CHECK and DELETE doubles)
prestim  = 500;
poststim = prestim + [1500 4000];  % cut ERPA shorter in 'phase' condition, in case of going for polynomial fit
pertwin = 5*srate; % total time window for ERPA
oswin = prestim + 3.3*srate;
rtwin = prestim + .8*srate;
% Define segment to analyze
segment = [prestim:poststim(1)];
ossegment = [oswin:length(avgBehavERI_mins{2})]; % over-shooting segment (post perturbation 'Stop', for tempo-change only)
rtsegment = [prestim:rtwin];
% Initialize integrals
[intAudio_mins,intAudio_null,intAudio_plus,intMotor_mins,intMotor_null,intMotor_plus,intBehav_mins,intBehav_null,intBehav_plus, ...
    overshoot_mins, overshoot_null, overshoot_plus, rt_mins, rt_plus]  =  deal(nan(nsubs,1)); % included overshooting post-tempo-change, and reaction time to phase-shift

% Initialize for plotting reaction times onsets
[pointer_mins,pointer_plus] = deal( nan( length(rtsegment), nsubs)  );
% Compute
figure(1), clf, hold on % to visualize fits
for subi = 1:nsubs
    
    % EEG (phase-shifts)
    intAudio_mins(subi) = -trapz(avgAudioERI_mins{1}(subi,segment)); % Note: to avoid trivial results on post-hoc comparisons, flip sign
    intAudio_null(subi) = trapz(avgAudioERI_null{1}(subi,segment));
    intAudio_plus(subi) = trapz(avgAudioERI_plus{1}(subi,segment));
    intMotor_mins(subi) = -trapz(avgMotorERI_mins{1}(subi,segment));
    intMotor_null(subi) = trapz(avgMotorERI_null{1}(subi,segment));
    intMotor_plus(subi) = trapz(avgMotorERI_plus{1}(subi,segment));
    
    % Behaviour (phase-shifts)
    intBehav_mins(subi) = -trapz(avgBehavERI_mins{1}(segment,subi));
    intBehav_null(subi) = trapz(avgBehavERI_null{1}(segment,subi));
    intBehav_plus(subi) = trapz(avgBehavERI_plus{1}(segment,subi));
    % Overshooting (tempo-change)
    overshoot_mins(subi) = -trapz(avgBehavERI_mins{2}(ossegment,subi));
    overshoot_plus(subi) = trapz(avgBehavERI_plus{2}(ossegment,subi));
    % Reaction times (phase)
    [tempSig_mins,~] = sigm_fit(rtsegment, avgBehavERI_mins{1}(rtsegment,subi), [],[],1); % you get inflection poin in ms, from perturbation onset
    if subi ~= 16 % hack: the model does not converge for sub#16... but it does if I flip sign; the inflexion point is invariant to this
    [tempSig_plus,~] = sigm_fit(rtsegment, avgBehavERI_plus{1}(rtsegment,subi), [],[],1); 
    else
    [tempSig_plus,~] = sigm_fit(rtsegment, -avgBehavERI_plus{1}(rtsegment,subi), [],[],1);
    end
        
    rt_mins(subi) = round( tempSig_mins(3) ) - prestim; % assign 3rd parameter (i.e., inflection point); round cause you need indexes on ms units
    rt_plus(subi) = round( tempSig_plus(3) ) - prestim; 
    % for plotting RTs
    pointer_mins(rt_mins(subi),subi) = avgBehavERI_mins{1}(rt_mins(subi),subi);
    pointer_plus(rt_plus(subi),subi) = avgBehavERI_plus{1}(rt_plus(subi),subi);
end


clearvars tempSig_mins tempSig_plus

% % Visualize RTs based on sigmoid fit
% for subi = 1:nsubs
%     figure(2); clf
% 
% hold on
% plot(avgBehavERI_mins{1}(rtsegment,subi), 'r' , 'LineWidth' , 1.5)
% plot(avgBehavERI_plus{1}(rtsegment,subi), 'g' , 'LineWidth' , 1.5)
% xline(rt_mins(subi),'r--', 'LineWidth' , 1.3)
% xline(rt_plus(subi),'g--', 'LineWidth' , 1.3)
% title(['Subject #' num2str(subi)])
% pause()
% end

% Assign outputs to R

% NOTE: by turning the '-' responses to negative, we control for eventual
% 'trivial results' in the model

eeg2r_phase = [intAudio_mins; intAudio_null; intAudio_plus; intMotor_mins; intMotor_null; intMotor_plus];

eeg2r_tempo = [prestim:poststim(2); ...
               -avgAudioERI_mins{2}(:,prestim:poststim(2)); avgAudioERI_null{2}(:,prestim:poststim(2)) ; avgAudioERI_plus{2}(:,prestim:poststim(2)); ... 
               -avgMotorERI_mins{2}(:,prestim:poststim(2)); avgMotorERI_null{2}(:,prestim:poststim(2)) ; avgMotorERI_plus{2}(:,prestim:poststim(2))];

bhv2r_phase = [intBehav_mins; intBehav_null; intBehav_plus];

bhv2r_tempo = [prestim:poststim(2); ...
               -avgBehavERI_mins{2}(prestim:poststim(2),:)'; avgBehavERI_null{2}(prestim:poststim(2),:)' ; avgBehavERI_plus{2}(prestim:poststim(2),:)'];
           
over2r_tempo = [overshoot_mins ; overshoot_plus];

rt2r_phase = [rt_mins ; rt_plus];
           
% Downsample output (NOTE: verify that regular intervals in timestamp indexes is not problematic for orthogonal polynomials)
dwn = 10; % downsampling factor
eeg2r_tempo = downsample(eeg2r_tempo(:,1:end-1)',dwn)';
bhv2r_tempo = downsample(bhv2r_tempo(:,1:end-1)',dwn)';
 

% Identify outliers     % TODO!: remove sub 6 from the zscores count (without removing, for not altering the indexes)
thresh = 2.23;

for condi = 1:2
    
    % EEG
    out_eegAudio_mins = zscore( mean(avgAudioERI_mins{condi},2) ); % note: all values should be signed (e.g: the biggest outliers would go in opposite +/- direction)
    out_eegAudio_null = zscore( mean(avgAudioERI_null{condi},2) );
    out_eegAudio_plus = zscore( mean(avgAudioERI_plus{condi},2) );
    out_eegMotor_mins = zscore( mean(avgMotorERI_mins{condi},2) );
    out_eegMotor_null = zscore( mean(avgMotorERI_null{condi},2) );
    out_eegMotor_plus = zscore( mean(avgMotorERI_plus{condi},2) );
    %Behaviour
    out_bhv_mins = zscore( mean(avgMotorERI_mins{condi},2) );
    out_bhv_null = zscore( mean(avgMotorERI_null{condi},2) );
    out_bhv_plus = zscore( mean(avgMotorERI_plus{condi},2) );
    %Reaction times
    out_rt_mins = zscore(rt_mins);
    out_rt_plus = zscore(rt_plus);

    % Display outliers
    disp([ condlabels{condi+5} ' Audio - : ' num2str( [find(abs(out_eegAudio_mins) > thresh)]' ) ]);
    disp([ condlabels{condi+5} ' Audio 0 : ' num2str( [find(abs(out_eegAudio_null) > thresh)]' ) ]);
    disp([ condlabels{condi+5} ' Audio + : ' num2str( [find(abs(out_eegAudio_plus) > thresh)]' ) ]);
    disp([ condlabels{condi+5} ' Motor - : ' num2str( [find(abs(out_eegMotor_mins) > thresh)]' ) ]);
    disp([ condlabels{condi+5} ' Motor 0 : ' num2str( [find(abs(out_eegMotor_null) > thresh)]' ) ]);
    disp([ condlabels{condi+5} ' Motor + : ' num2str( [find(abs(out_eegMotor_plus) > thresh)]' ) ]);
    disp([ condlabels{condi+5} ' Behaviour - : ' num2str( [find(abs(out_bhv_mins) > thresh)]' ) ]);
    disp([ condlabels{condi+5} ' Behaviour 0 : ' num2str( [find(abs(out_bhv_null) > thresh)]' ) ]);
    disp([ condlabels{condi+5} ' Behaviour + : ' num2str( [find(abs(out_bhv_plus) > thresh)]' ) ]);
    disp([ condlabels{condi+5} ' RT - : ' num2str( [find(abs(out_rt_mins) > thresh)]' ) ]);
    disp([ condlabels{condi+5} ' RT + : ' num2str( [find(abs(out_rt_plus) > thresh)]' ) ]);
    
end

% Export to R
cd('/Users/mattiaipem/Documents/IPEM/Projects/ Lou _ Entrainment/ERFA_Metronomes/Data/Experiment');           
writematrix(eeg2r_phase , 'temp_eegERFA_Phase.xlsx')
writematrix(eeg2r_tempo , 'temp_eegERFA_Tempo.xlsx') 
writematrix(bhv2r_phase , 'temp_bhvERFA_Phase.xlsx')
writematrix(bhv2r_tempo , 'temp_bhvERFA_Tempo.xlsx') 
writematrix(over2r_tempo , 'temp_overshoot_Tempo.xlsx') 
writematrix(rt2r_phase , 'temp_rt_Phase.xlsx') 
