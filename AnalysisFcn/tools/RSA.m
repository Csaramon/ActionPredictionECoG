%(1) Epoch the data around your timepoint of interest:


eeg_epoched                   = []; % get data into fieldtrip format

eeg_epoched.fsample    = eeg_resampled.fsample; % sampling frequency

eeg_epoched.label          = eeg_resampled.label; % channel labels

% loop through trials (for each trial, you need a channel x timepoints
% matrix)
for n_trial = 1:size(beh.cue.allsampleidx, 1)
        
    % time region of interest (including a baseline period)
    twoi    = [beh.iti_last2sec.allsampleidx{n_trial, 1}, beh.cue.allsampleidx{n_trial, 1}];
    fprintf('Epoch trial %d: %d --> %d (%d samples)...\n', n_trial, twoi(1), twoi(end), size(twoi, 2));
        
    eeg_epoched.trial{n_trial, 1}  = eeg_resampled.trial{1}(:, twoi); % electrodes x timepoints within one trial
    eeg_epoched.time{n_trial, 1}   = -2 : 1/eeg_resampled.fsample : ((size(twoi, 2) - 1)/eeg_resampled.fsample - 2); % in [sec]
end

%(2) Perform ICA to transform the EEG data into information sources (using fieldtrip)


cfg                 = [];
cfg.method          = 'runica';
cfg.channel     = subjectdata.usableelectrodes;
cfg.trials          = 'all'; % or a selection given as a 1xN vector (default = 'all')
cfg.numcomponent    = 'all'; % or number (default = 'all')
cfg.demean          = 'yes'; % or 'no', whether to demean the input data (default = 'yes')
cfg.updatesens      = 'yes'; % or 'no' (default = 'yes')
cfg.randomseed      = 1;
        
% perform ICA
eeg_components      = ft_componentanalysis(cfg, eeg_epoched);
%(3) Perform baseline correction in each trial (using fieldtrip)


cfg                 = [];
cfg.continuous      = 'no';
cfg.channel         = 'all';
cfg.demean          = 'yes'; % whether to perform baseline correction
cfg.baselinewindow  = [-.2, 0]; % time window for baseline correction
    
eeg_baseline        = ft_preprocessing(cfg, eeg_epoched);

%(4) Select analysis-relevant timepoints


cfg         = [];
cfg.trials  = 'all';
cfg.toilim  = [-.5, 2]; % to specify a latency window in seconds
    
eeg_twoi    = ft_redefinetrial(cfg, eeg_baseline);

%(5) Re-arrange data for easier access during the RSA


% get data into correct format
all_data = nan(size(eeg_twoi.trial, 2), length(eeg_twoi.label), ...
    size(eeg_twoi.time{1}, 2)); % trials x components x timepoints-within-trial
for n_trial = 1:size(eeg_twoi.trial, 2)
    all_data(n_trial, :, :)     = eeg_twoi.trial{n_trial};
end

objects                     = beh.objects; % one specific object per trial, this is the relevant behavioral information for the RSA


%(6) Perform the RSA


% preallocate interesting output
on_diagonal_RSA                     = nan(size(all_data, 3), 1); % on-diagonal-RSA-values
off_diagonal_RSA                    = nan(size(all_data, 3), 1); % off-diagonal-RSA-values
    
% loop through timepoints
for n_timepoint = 1:size(all_data, 3)
        
    % give info how far you are
    if mod(n_timepoint, 500) == 0
        fprintf('Analysing timepoint %d ...\n', n_timepoint);
    end
        
    if b_identicalRandAllTp % e.g., b_identicalRandAllTp = true
        rng(isub); % isub is the subject index
    end
        
    % select ERP-data across channels for this timepoint (temporal averaging)
    timepoints2select   = n_timepoint - temporalsmooth : n_timepoint + temporalsmooth; % e.g., temporalsmooth = 15 at a sampling frequency of 1000 Hz
    timepoints2select   = timepoints2select(timepoints2select > 0 & timepoints2select <= size(all_data, 3)); % correction at the borders of the analysis time window
    this_data           = mean(all_data(:, :, timepoints2select), 3); % trials x features
        
    % preallocate the confusion matrix (pairwise correlations between
    % all object representations)
    tmp_confusion_matrix    = nan(size(unique(objects), 1), size(unique(objects), 1), numrepetitions); % nobjects, nobjects, repetitions
        
    % perform several repetitions of the analysis to make it
    % independent from a specific half-selection
    for n_rep = 1:numrepetitions % e.g., numrepetitions = 10
            
        % distribute onto different halves
        b_allobjectsinbothhalves = false;
        while b_allobjectsinbothhalves == false
                
            % distribute overall data randomly onto both data halves
            randtmp         = rand(size(this_data, 1), 1);
            halfidx         = randtmp > median(randtmp);
            this_data_half1 = this_data(halfidx == 1, :);
            objects_half1   = objects(halfidx == 1, :);
            this_data_half2 = this_data(halfidx == 0, :);
            objects_half2   = objects(halfidx == 0, :);
                
            % check whether all objects were sampled in both data halves
            if length(unique(objects_half1)) == length(unique(objects)) && length(unique(objects_half2)) == length(unique(objects))
                b_allobjectsinbothhalves = true;
            end
        end
            
        sel_data_half1  = this_data_half1; % ntrials x features
        sel_data_half2  = this_data_half2;
            
        % preallocate
        representation_half1 = nan(size(sel_data_half1, 2), size(unique(objects), 1)); % features x nobjects
        representation_half2 = nan(size(representation_half1));
            
        % loop through unique object numbers
        unique_objects  = unique(objects);
        for n_obj = 1:size(unique_objects, 1)
                
            % get "name" of this object (0:7)
            obj_number      = unique_objects(n_obj); % convert object-index to object-name
                
            % average ERPs to obtain representation (across trials)
            representation_half1(:, n_obj) = mean(sel_data_half1(objects_half1 == obj_number, :), 1);
            representation_half2(:, n_obj) = mean(sel_data_half2(objects_half2 == obj_number, :), 1);
        end
            
        % calculate confusion matrix and perform
        % transform correlation values using Fisher-z-transformation
        tmp_confusion_matrix(:, :, n_rep)   = atanh(corr(representation_half1, representation_half2, 'Type', corr_type)); % e.g., corr_type = 'Spearman'
       
    end % end of repetition-loop
        
    % average confusion matrix (make the result independent from a
    % specific assignment of trials to halves)
    confusion_matrix    = mean(tmp_confusion_matrix, 3); % average across repetitions
   
    % extract RSA values for identical items
    on_diagonal_fields                  = logical(eye(size(confusion_matrix)));
    on_diagonal_RSA(n_timepoint, 1)     = mean(confusion_matrix(on_diagonal_fields)); % similarity of identical objects
           
    % extract RSA values for non-identical items
    tmp_idx                 = randperm(size(confusion_matrix, 1));
    normal_idx              = 1:size(unique_objects, 1);
    while any((tmp_idx' - normal_idx') == 0)
        tmp_idx             = randperm(size(confusion_matrix, 1)); % ensure that no on-diagonal-field is contained
    end
    off_diagonal_fields                 = on_diagonal_fields(tmp_idx, :);
    off_diagonal_RSA(n_timepoint, 1)    = mean(confusion_matrix(off_diagonal_fields)); % similarity of different objects
end % end of timepoint-loop
