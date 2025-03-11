function [ conmat ] = rand_SSVEP_workmem(p,RDK,flag_training)
%rand_FShiftBase randomizes experimental conditions
% move onset only works for constant frequency for all RDKs (i.e. 120)

% open issues:
% if we want to randomize the post-cue tracking time window (which we have to do to make sure that participants are tracking
% constantly throughout the trial), trials with the current frequencies are only allowed by trial-onset and they don't sync
% throughout the trial. Hence, it's only possible to analyze SSVEPs in a post-cue time window (i.e. from 500 to 1500 ms) but not
% in a pre-retention time-window (i.e. from -1000 to 0) as these signals are not phase locked


% set trial number etc
if flag_training~=0
    conmat.totaltrials = p.stim.trialnum_train;
    conmat.totalblocks = 1;
else
    conmat.totaltrials =  sum(p.stim.con_repeats+p.stim.con_repeats_catch);
    conmat.totalblocks = p.stim.blocknum;
end
conmat.trialsperblock = conmat.totaltrials/conmat.totalblocks;

% matrix with onset times of on framesfor RDKs
[t.onframesonset, t.onframesonset_times]= deal(nan(numel(RDK.RDK),p.scr_refrate*max(p.stim.time_postcuetrack)));
for i_rdk = 1:numel(RDK.RDK)
    t.mat = ceil(1:p.scr_refrate/RDK.RDK(i_rdk).freq:size(t.onframesonset,2));
    t.onframesonset(i_rdk,t.mat)=1;
    t.onframesonset_times(i_rdk,t.mat)=t.mat./p.scr_refrate;
end

% move
[t.movonset_frames, t.movonset_times]=deal(nan(1,p.scr_refrate*max(p.stim.time_postcuetrack)));
t.mat = 1:p.scr_refrate/RDK.RDK(1).mov_freq:size(t.movonset_frames,2);
t.movonset_frames(t.mat)=1;
t.movonset_times(t.mat)=t.mat./p.scr_refrate;


%% start randomization
% randomize catch or regular trials [1=regular, 2=catch]
if flag_training~=0 % training
    conmat.mats.trialtype = randsample([repmat(1,1,sum(p.stim.con_repeats)) repmat(2,1,sum(p.stim.con_repeats_catch))],conmat.totaltrials);
else
    conmat.mats.trialtype = [repmat(1,1,sum(p.stim.con_repeats)) repmat(2,1,sum(p.stim.con_repeats_catch))];
end

% randomize condition
if flag_training~=0 % training
    t.mat = [repmat(p.stim.condition,1,floor(conmat.totaltrials/numel(p.stim.condition))) ...
        randsample(p.stim.condition,mod(conmat.totaltrials, numel(p.stim.condition)))];
    conmat.mats.condition = t.mat(randsample(numel(t.mat),numel(t.mat)));
else
    conmat.mats.condition = [ ... % for regular and catch trials create con labels according to repetition parameter
        cell2mat(arrayfun(@(x,y) repmat(x,1,y), p.stim.condition, p.stim.con_repeats, 'UniformOutput', false)), ...
        cell2mat(arrayfun(@(x,y) repmat(x,1,y), p.stim.condition, p.stim.con_repeats_catch, 'UniformOutput', false))];
end


% assign cue (1 = attend to RDK 1 (left); 2 = attend to RDK 4 (right))
conmat.mats.cue = nan(1,conmat.totaltrials);
conmat.mats.cue(ismember(conmat.mats.condition,1:2:6))=1;
conmat.mats.cue(ismember(conmat.mats.condition,2:2:6))=2;


% randomize eventtype (1 = probe match; 2 = probe mismatch)
conmat.mats.eventtype = nan(1,conmat.totaltrials);
t.evtype = [repmat(1,1,p.stim.event.ratio(1)) repmat(2,1,p.stim.event.ratio(2))]; % pre-assign values according to ratio of event types
for i_trtype = 1:2 % loop across catch and regular trials
    t.idx1 = conmat.mats.trialtype == i_trtype;
    for i_con = 1:numel(p.stim.condition)
        t.idx2 = conmat.mats.condition==p.stim.condition(i_con);
        for i_cue = unique(conmat.mats.cue)
            t.idx3 = conmat.mats.cue == i_cue;
            % distribute equally across condition and cue to left or right
            t.evnum = sum(t.idx1&t.idx2&t.idx3);
            % fill up as many slots as possible uniformly, random sample for the rest
            t.evtypemat = [repmat(t.evtype, 1, floor(t.evnum/numel(t.evtype))) randsample(t.evtype, mod(t.evnum,numel(t.evtype)))];

            conmat.mats.eventtype(t.idx1&t.idx2&t.idx3) = t.evtypemat;
        end
    end
end


% determine event RDK
conmat.mats.eventRDK = nan(1,conmat.totaltrials);
% for con 1,2,3,4 left or right RDK1 or RDK4
conmat.mats.eventRDK(ismember(conmat.mats.condition,[1 5])) = 1;
conmat.mats.eventRDK(ismember(conmat.mats.condition,[2 6])) = 4;
% for conditions 3 and 4, event change could be on RDK1,2 or RDK 4 and RDK 5
for i_trtype = 1:2 % loop across catch and regular trials
    for i_evtype = 1:2
        % for con 3 RDK
        t.idx = conmat.mats.trialtype==i_trtype & conmat.mats.eventtype==i_evtype & conmat.mats.condition==3;
        conmat.mats.eventRDK(t.idx) = [repmat([1 2],1,floor(sum(t.idx)/2)) randsample([1 2], mod(sum(t.idx), 2))];
        t.idx = conmat.mats.trialtype==i_trtype & conmat.mats.eventtype==i_evtype & conmat.mats.condition==4;
        conmat.mats.eventRDK(t.idx) = [repmat([4 5],1,floor(sum(t.idx)/2)) randsample([4 5], mod(sum(t.idx), 2))];
    end
end

% randomize orientations of stimuli
conmat.mats.stimangles = nan([numel(RDK.RDK),RDK.RDK(1).num,conmat.totaltrials]);
conmat.mats.stimangles(:) = randsample(p.stim.angles,numel(conmat.mats.stimangles),true);

% randomize orientation rotations of stimuli
conmat.mats.eventangles = conmat.mats.stimangles;
t.idx1 = conmat.mats.eventtype==2; % only for non match stimuli, one stimulus will be changed
t.rdks = unique(conmat.mats.eventRDK); % which RDk will contain event?
for i_rdk = t.rdks
    % index event for certain rdk and change as event type
    t.idx2 = conmat.mats.eventRDK==i_rdk & t.idx1;
    % randomly sample one of the event stimuli
    t.eventstim = randsample(RDK.RDK(i_rdk).num,sum(t.idx2),true);
    t.eventrotation = randsample(p.stim.event.rotations,sum(t.idx2),true);
    % bring it all together
    t.idxmat = false(size(conmat.mats.eventangles));
    % sub2ind finds the subindices defined by each vector
    t.idxmat(sub2ind(size(t.idxmat), repmat(i_rdk,1,sum(t.idx2)), t.eventstim', find(t.idx2))) = true;
    
    % now change values by rotation values
    conmat.mats.eventangles(t.idxmat) = conmat.mats.eventangles(t.idxmat)+t.eventrotation';
end

% randomize post-cue times [across trialtype and condition]
t.times = {p.stim.time_postcuetrack(1):0.1:p.stim.time_postcuetrack(2);
    p.stim.time_postcuetrack_catch(1):0.1:p.stim.time_postcuetrack_catch(2)-0.01};
conmat.mats.postcue_times = nan(1,conmat.totaltrials);

for i_trtype = 1:2 % loop across catch and regular trials
    t.idx1 = conmat.mats.trialtype == i_trtype; % for regular trials
    for i_con = 1:numel(p.stim.condition) % loop across conditions
        t.idx2 = conmat.mats.condition==p.stim.condition(i_con);
        % uniformly distribute all possible times and randsample the rest
        t.mat = [repmat(t.times{i_trtype},1,floor(sum(t.idx1&t.idx2)/numel(t.times{i_trtype}))), ...
            randsample(t.times{i_trtype}, mod(sum(t.idx1&t.idx2), numel(t.times{i_trtype})))];
        conmat.mats.postcue_times(t.idx1&t.idx2) = t.mat(randsample(sum(t.idx1&t.idx2),sum(t.idx1&t.idx2)));
    end
end

conmat.mats.postcue_frames = round(conmat.mats.postcue_times*p.scr_refrate);

% randomize retention times [across trialtype and condition]
t.times = p.stim.time_retention(1):0.1:p.stim.time_retention(2);
conmat.mats.retention_times = nan(1,conmat.totaltrials);

for i_trtype = 1:2 % loop across catch and regular trials
    t.idx1 = conmat.mats.trialtype == i_trtype; % for regular trials
    for i_con = 1:numel(p.stim.condition) % loop across conditions
        t.idx2 = conmat.mats.condition==p.stim.condition(i_con);
        % uniformly distribute all possible times and randsample the rest
        t.mat = [repmat(t.times,1,floor(sum(t.idx1&t.idx2)/numel(t.times))), ...
            randsample(t.times, mod(sum(t.idx1&t.idx2), numel(t.times)))];
        conmat.mats.retention_times(t.idx1&t.idx2) = t.mat(randsample(sum(t.idx1&t.idx2),sum(t.idx1&t.idx2)));
    end
end
conmat.mats.retention_frames = round(conmat.mats.retention_times*p.scr_refrate);


% randomize pre-cue times
t.times = p.stim.time_precue(1):0.1:p.stim.time_precue(2);
conmat.mats.precue_times = nan(1,conmat.totaltrials);

for i_trtype = 1:2 % loop across catch and regular trials
    t.idx1 = conmat.mats.trialtype == i_trtype; % for regular trials
    for i_con = 1:numel(p.stim.condition) % loop across conditions
        t.idx2 = conmat.mats.condition==p.stim.condition(i_con);
        % uniformly distribute all possible times and randsample the rest
        t.mat = [repmat(t.times,1,floor(sum(t.idx1&t.idx2)/numel(t.times))), ...
            randsample(t.times, mod(sum(t.idx1&t.idx2), numel(t.times)))];
        conmat.mats.precue_times(t.idx1&t.idx2) = t.mat(randsample(sum(t.idx1&t.idx2),sum(t.idx1&t.idx2)));
    end
end

conmat.mats.precue_frames = round(conmat.mats.precue_times*p.scr_refrate);

% include pre-cue events
% number of pre-cue events
t.mat = [repmat(p.stim.precue_event.num,1,floor(conmat.totaltrials/numel(p.stim.precue_event.num))), ...
    randsample(p.stim.precue_event.num,mod(conmat.totaltrials, numel(p.stim.precue_event.num)))];
conmat.mats.precue_eventnum = randsample(t.mat,numel(t.mat));

% indicate event: 1 - left shorter, 2 - left longer, 3 - right shorter, 4 - right longer
conmat.mats.precue_eventid = nan(size(conmat.mats.precue_eventnum));
t.mat = [repmat([1 2 3 4],1,floor(sum(conmat.mats.precue_eventnum)/4)),randsample([1 2 3 4],mod(sum(conmat.mats.precue_eventnum), 4))];
conmat.mats.precue_eventid(conmat.mats.precue_eventnum==1)= t.mat(randsample(numel(t.mat ),numel(t.mat )));

% target or distractor? 1 = target, 2= distractor
conmat.mats.precue_eventtype = nan(size(conmat.mats.precue_eventnum));
conmat.mats.precue_eventtype(ismember(conmat.mats.precue_eventid,p.stim.precue_event.targets))=1;
conmat.mats.precue_eventtype(ismember(conmat.mats.precue_eventid,setdiff([1 2 3 4], p.stim.precue_event.targets)))=2;

% pre-cue event onset times
t.idx = find(conmat.mats.precue_eventnum==1);
conmat.mats.precue_event_onset_times = nan(size(conmat.mats.precue_eventnum));
for i_tr = 1:numel(t.idx)
    % index possible positions
    t.idx1 = p.stim.precue_event.min_onset:0.1:conmat.mats.precue_times(t.idx(i_tr))-(p.stim.precue_event.length+p.stim.precue_event.min_offset);
    conmat.mats.precue_event_onset_times(t.idx(i_tr)) = randsample(t.idx1,1);
end
conmat.mats.precue_event_onset_frames = round(conmat.mats.precue_event_onset_times*p.scr_refrate);

%% randomize all information across experiment
t.tidx = randperm(conmat.totaltrials);
conmat.mats.trialtype = conmat.mats.trialtype(:,t.tidx);
conmat.mats.condition = conmat.mats.condition(:,t.tidx);
conmat.mats.cue = conmat.mats.cue(:,t.tidx);
conmat.mats.eventtype = conmat.mats.eventtype(:,t.tidx);
conmat.mats.eventRDK = conmat.mats.eventRDK(:,t.tidx);
conmat.mats.stimangles = conmat.mats.stimangles(:,:,t.tidx);
conmat.mats.eventangles = conmat.mats.eventangles(:,:,t.tidx);
conmat.mats.postcue_frames = conmat.mats.postcue_frames(:,t.tidx);
conmat.mats.postcue_times = conmat.mats.postcue_times(:,t.tidx);
conmat.mats.retention_frames = conmat.mats.retention_frames(:,t.tidx);
conmat.mats.retention_times = conmat.mats.retention_times(:,t.tidx);
conmat.mats.precue_frames = conmat.mats.precue_frames(:,t.tidx);
conmat.mats.precue_times = conmat.mats.precue_times(:,t.tidx);
conmat.mats.precue_eventnum = conmat.mats.precue_eventnum(:,t.tidx);
conmat.mats.precue_eventid = conmat.mats.precue_eventid(:,t.tidx);
conmat.mats.precue_eventtype = conmat.mats.precue_eventtype(:,t.tidx);
conmat.mats.precue_event_onset_times = conmat.mats.precue_event_onset_times(:,t.tidx);
conmat.mats.precue_event_onset_frames = conmat.mats.precue_event_onset_frames(:,t.tidx);

conmat.mats.block = repmat(1:conmat.totalblocks,conmat.trialsperblock,1);
conmat.mats.block = conmat.mats.block(:)';

%% write all information into trial structure
% create frame mat, onset time for events

for i_tr = 1:conmat.totaltrials
    % trialnumber
    conmat.trials(i_tr).trialnum = i_tr;

    % block number
    conmat.trials(i_tr).blocknum = conmat.mats.block(i_tr);

    % trial type
    conmat.trials(i_tr).trialtype = conmat.mats.trialtype(i_tr);
    
    % condition
    % [cue_left_3t3d cue_right_3t3d cue_left_6t cue_right_6t cue_left_3t cue_right_3t]
    % [RDK1 RDK3; RDK4 RDK6] [RDK1 RDK3; RDK4 RDK6] [RDK1 RDK2; RDK4 RDK5] [RDK1 RDK2; RDK4 RDK5] [RDK1; RDK4] [RDK1; RDK4] 
    conmat.trials(i_tr).condition = conmat.mats.condition(i_tr);
    
    % RDK to display
    t.mat = {[1 3;4 6]; [1 3;4 6]; [1 2; 4 5]; [1 2; 4 5]; [1; 4]; [1; 4]};
    conmat.trials(i_tr).RDK2display = t.mat{conmat.mats.condition(i_tr)};
    
    % cue ((RDK1, RDK2) [1,2])
    conmat.trials(i_tr).cue = conmat.mats.cue(i_tr);
    
    % type of events ((match, nonmatch) [1, 2])
    conmat.trials(i_tr).eventtype = conmat.mats.eventtype(:,i_tr);
    
    % which RDK shows event?
    conmat.trials(i_tr).eventRDK = conmat.mats.eventRDK(:,i_tr);
    
    % stim angles of the displayed stimuli
    conmat.trials(i_tr).stimangles =  conmat.mats.stimangles(:,:,i_tr);
    
    % event angles of the displayed stimuli
    conmat.trials(i_tr).eventangles =  conmat.mats.eventangles(:,:,i_tr);
    
    % pre-cue frames
    conmat.trials(i_tr).precue_frames = conmat.mats.precue_frames(:,i_tr);
    
    % pre-cue times
    conmat.trials(i_tr).precue_times = conmat.mats.precue_times(:,i_tr);
    
    % post-cue times
    conmat.trials(i_tr).postcue_times = conmat.mats.postcue_times(:,i_tr);
    
    % post-cue frames
    conmat.trials(i_tr).postcue_frames = conmat.mats.postcue_frames(:,i_tr);

    % retention times
    conmat.trials(i_tr).retention_times = conmat.mats.retention_times(:,i_tr);
    
    % retention frames
    conmat.trials(i_tr).retention_frames = conmat.mats.retention_frames(:,i_tr);

    % number of pre-cue events
    conmat.trials(i_tr).precue_eventnum = conmat.mats.precue_eventnum(:,i_tr);

    % what type of event [1 - left shorter, 2 - left longer, 3 - right shorter, 4 - right longer]
    conmat.trials(i_tr).precue_eventid = conmat.mats.precue_eventid(:,i_tr);

    % what type of event [1 - target; 2 - distractor]
    conmat.trials(i_tr).precue_eventtype = conmat.mats.precue_eventtype(:,i_tr);

    % times of pre-cue event
    conmat.trials(i_tr).precue_event_onset_times = conmat.mats.precue_event_onset_times(:,i_tr);

    % frames of pre-cue event
    conmat.trials(i_tr).precue_event_onset_frames = conmat.mats.precue_event_onset_frames(:,i_tr);
end

%% graphical checks
% % how many regular trials
% cell2mat(arrayfun(@(x) [sum(conmat.mats.condition==x & conmat.mats.trialtype==1); sum(conmat.mats.condition==x & conmat.mats.trialtype==2)], ...
%     unique(conmat.mats.condition), 'UniformOutput', false))
% 
% % events in all conditions distributed equally
% cell2mat(arrayfun(@(x) [sum(conmat.mats.condition==x & conmat.mats.trialtype==1); sum(conmat.mats.eventtype(conmat.mats.condition==x & conmat.mats.trialtype==1)==1)], ...
%     unique(conmat.mats.condition), 'UniformOutput', false))
% 
% % which RDK across conditions
% figure;
% for i_con = unique(conmat.mats.condition)
%     subplot(1, numel(unique(conmat.mats.condition)),i_con)
%     t.dat = conmat.mats.eventRDK(conmat.mats.condition==i_con & conmat.mats.trialtype==1);
%     histogram(t.dat,1:6)
%     title(sprintf("con %1.0f",i_con))
% end
% 
% % histogram of angles
% figure; histogram(conmat.mats.stimangles(:))
% 
% % histogram of stim angle modulations
% t.dat = conmat.mats.stimangles(:)-conmat.mats.eventangles(:); t.dat = t.dat(t.dat~=0);
% figure; histogram(t.dat)
% 
% % histogram of post-cue retention times
% figure; histogram(conmat.mats.postcue_times(:),numel(unique(conmat.mats.postcue_times)))
% 
% % histogram of pre-cue times
% figure; histogram(conmat.mats.precue_times(:),numel(unique(conmat.mats.precue_times)))
% 
% % histogram of pre-cue times
% figure; histogram(conmat.mats.retention_times(:),numel(unique(conmat.mats.retention_times)))
% 
% % total stim time
% sum(conmat.mats.precue_times+conmat.mats.postcue_times+conmat.mats.retention_times+1+1.5)/60

    

end

