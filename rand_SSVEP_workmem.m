function [ conmat ] = rand_SSVEP_workmem(p,RDK,flag_training)
%rand_FShiftBase randomizes experimental conditions
% move onset only works for constant frequency for all RDKs (i.e. 120)




% set trial number etc
if flag_training~=0
    conmat.totaltrials = sum(p.stim.con_repeats+p.stim.con_repeats_catch);
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
% randomize condition
t.mat = repmat(p.stim.condition,conmat.totaltrials/numel(p.stim.condition),1);
conmat.mats.condition = t.mat(:)';


% assign cue (1 = attend to RDK 1 (left); 2 = attend to RDK 4 (right))
conmat.mats.cue = nan(1,conmat.totaltrials);
conmat.mats.cue(ismember(conmat.mats.condition,1:2:6))=1;
conmat.mats.cue(ismember(conmat.mats.condition,2:2:6))=2;


% randomize eventtype (1 = probe match; 2 = probe mismatch)
conmat.mats.eventtype = nan(1,conmat.totaltrials);
t.evtype = [repmat(1,1,p.stim.event.ratio(1)) repmat(2,1,p.stim.event.ratio(2))]; % pre-assign values according to ratio of event types
for i_con = 1:numel(p.stim.condition)
    t.idx = conmat.mats.condition==p.stim.condition(i_con);
    for i_cue = 1:2
        t.idx2 = conmat.mats.cue == i_cue;
        % distribute equally across condition and cue to left or right
        t.evnum = sum(t.idx&t.idx2);
        % fill up as many slots as possible uniformly, random sample for the rest
        t.evtypemat = [repmat(t.evtype, 1, floor(t.evnum/numel(t.evtype))) randsample(t.evtype, mod(t.evnum,numel(t.evtype)))];

        conmat.mats.eventtype(t.idx&t.idx2) = t.evtypemat;
    end
end


% determine event RDK
conmat.mats.eventRDK = nan(1,conmat.totaltrials);
% for con 1,2,3,4 left or right RDK1 or RDK4
conmat.mats.eventRDK(ismember(conmat.mats.condition,[1 5])) = 1;
conmat.mats.eventRDK(ismember(conmat.mats.condition,[2 6])) = 4;
% for conditions 3 and 4, event change could be on RDK1,2 or RDK 4 and RDK 5
for i_evtype = 1:2
    % for con 3 RDK
    t.idx = conmat.mats.eventtype==i_evtype & conmat.mats.condition==3;
    conmat.mats.eventRDK(t.idx) = [repmat([1 2],1,floor(sum(t.idx)/2)) randsample([1 2], mod(sum(t.idx), 2))];
    t.idx = conmat.mats.eventtype==i_evtype & conmat.mats.condition==4;
    conmat.mats.eventRDK(t.idx) = [repmat([4 5],1,floor(sum(t.idx)/2)) randsample([4 5], mod(sum(t.idx), 2))];
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


%% include pre-cue events?
%% randomize pre-cue times
%% randomize post-cue times


%%
% pre-allocate possible presentation times
conmat.mats.event_onset_frames = nan(max(p.stim.eventnum),conmat.totaltrials);
t.poss_frames = find(p.stim.event.min_onset<t.movonset_times & ...
    t.movonset_times<(p.stim.time_postcue-p.stim.event.length-p.stim.event.min_offset));
t.poss_frames_1 = find(p.stim.event.min_onset<t.movonset_times & ...
    t.movonset_times<(p.stim.time_postcue-p.stim.event.length-p.stim.event.min_offset-p.stim.event.min_dist-0.01));
t.poss_frames_2 = find(p.stim.event.min_onset+p.stim.event.min_dist<t.movonset_times & ...
    t.movonset_times<(p.stim.time_postcue-p.stim.event.length-p.stim.event.min_offset));
% loop across conditions
for i_con = 1:numel(p.stim.condition)
    % for single events first
    t.idx = repmat(conmat.mats.cue,2,1) == i_con & repmat(conmat.mats.eventnum,2,1) == 1 & ~isnan(conmat.mats.eventtype);
    if sum(t.idx(:))<numel(t.poss_frames) % if more possible positons than actual events
        t.poss_frames_mat = t.poss_frames(randsample(numel(t.poss_frames),sum(t.idx(:))));
    else
        t.poss_frames_mat = [repmat(t.poss_frames,1,floor(sum(t.idx(:))/numel(t.poss_frames)))...
            t.poss_frames(randsample(numel(t.poss_frames),mod(sum(t.idx(:)),numel(t.poss_frames))))];
    end
    % write to frame mat
    conmat.mats.event_onset_frames(t.idx)=t.poss_frames_mat;
    
    % for two events
    t.idx = repmat(conmat.mats.cue,2,1) == i_con & repmat(conmat.mats.eventnum,2,1) == 2;
    t.idx2=find(t.idx(1,:));
    t.poss_frames_mat = [];
    % loop across all events
    for i_ev = 1:numel(t.idx2)
        t.poss_frames_mat(1,i_ev)= t.poss_frames_1(randperm(numel(t.poss_frames_1),1));
        % index all positions still possible
        t.idx3 = t.poss_frames_2>(t.poss_frames_mat(1,i_ev)+p.stim.event.min_dist*p.scr_refrate);
        t.idx4 = find(t.idx3);
        % randomly select from possible frames (uniqform distribution)
        t.idx5 = t.idx4(randsample(numel(t.idx4),1));
        % randomly select from pssoble frames (beta distribution random number) --> compensate for righward distribution?
%         t.idx5 = round(betarnd(1,5,1)*(t.idx4(end)-t.idx4(1))+t.idx4(1));
        t.idx5 = round(betarnd(1.2,3,1)*(t.idx4(end)-t.idx4(1))+t.idx4(1));
        t.poss_frames_mat(2,i_ev)=t.poss_frames_2(t.idx5);
    end
    conmat.mats.event_onset_frames(t.idx)=t.poss_frames_mat;
    
end
conmat.mats.event_onset_times = conmat.mats.event_onset_frames./p.scr_refrate;
% % graphical check
% figure; subplot(2,1,1);histogram(conmat.mats.event_onset_frames(:),50);subplot(2,1,2);histogram(conmat.mats.event_onset_times(:),50)
% figure; subplot(2,1,1);histogram(diff(conmat.mats.event_onset_frames),50);subplot(2,1,2);histogram(conmat.mats.event_onset_times(:),50)
% 
% for i_tr = 1:100
% test(i_tr,:,:) = conmat.mats.event_onset_times;
% end
% figure; subplot(2,1,1); histogram(test(:)); subplot(2,1,2); histogram(diff(test,1,2))


% randomize pre-cue times
% all possible pre_cue_frames
t.allframes = p.stim.time_precue(1)*p.scr_refrate:p.stim.time_precue(2)*p.scr_refrate;
t.allframes = t.allframes(mod(t.allframes,p.scr_imgmultipl)==0); % only frames that are integers of frames per flip (i.e. 4)
if conmat.totaltrials<numel(t.allframes)
    conmat.mats.pre_cue_frames = t.allframes(randsample(1:numel(t.allframes),conmat.totaltrials));
else
    conmat.mats.pre_cue_frames = [repmat(t.allframes,1,floor(conmat.totaltrials/numel(t.allframes))) ...
        t.allframes(round(linspace(1,numel(t.allframes),mod(conmat.totaltrials,numel(t.allframes)))))];
end
conmat.mats.pre_cue_frames = conmat.mats.pre_cue_frames(randperm(numel(conmat.mats.pre_cue_frames)));
conmat.mats.pre_cue_times = conmat.mats.pre_cue_frames./p.scr_refrate;

% add pre-cue frames to events
conmat.mats.event_onset_times = conmat.mats.event_onset_times+conmat.mats.pre_cue_times;
conmat.mats.event_onset_frames = conmat.mats.event_onset_frames + conmat.mats.pre_cue_frames;

%% randomize all information across experiment
t.tidx = randperm(conmat.totaltrials);
conmat.mats.condition = conmat.mats.condition(:,t.tidx);
conmat.mats.cue = conmat.mats.cue(:,t.tidx);
conmat.mats.eventnum = conmat.mats.eventnum(:,t.tidx);
conmat.mats.eventtype = conmat.mats.eventtype(:,t.tidx);
conmat.mats.eventRDK = conmat.mats.eventRDK(:,t.tidx);
conmat.mats.eventdirection = conmat.mats.eventdirection(:,t.tidx);
conmat.mats.event_onset_frames = conmat.mats.event_onset_frames(:,t.tidx);
conmat.mats.event_onset_times = conmat.mats.event_onset_times(:,t.tidx);
conmat.mats.pre_cue_frames = conmat.mats.pre_cue_frames(:,t.tidx);
conmat.mats.pre_cue_times = conmat.mats.pre_cue_times(:,t.tidx);

conmat.mats.block = repmat(1:conmat.totalblocks,conmat.trialsperblock,1);
conmat.mats.block = conmat.mats.block(:);

%% write all information into trial structure
% create frame mat, onset time for events

for i_tr = 1:conmat.totaltrials
    % trialnumber
    conmat.trials(i_tr).trialnum = i_tr;
    
    % block number
    conmat.trials(i_tr).blocknum = conmat.mats.block(i_tr);
    
    % condition [C1 C2; C1 C2] [C1 C2; C1 C2] [C1 C2; C1 C3] [C1 C2; C1 C3] [C1 C2; C2 C3] [C1 C2; C2 C3]
    conmat.trials(i_tr).condition = conmat.mats.condition(i_tr);
    
    % RDK to display
    t.mat = [p.stim.RDKcenter p.stim.RDKperi+2];
    conmat.trials(i_tr).RDK2display = t.mat(conmat.mats.condition(i_tr),:);
    
    % cue ((RDK1, RDK2) [1,2])
    conmat.trials(i_tr).cue = conmat.mats.cue(i_tr);
    
    % number of events
    conmat.trials(i_tr).eventnum = conmat.mats.eventnum(i_tr);
    
    % type of events ((target, distractor) [1, 2])
    conmat.trials(i_tr).eventtype = conmat.mats.eventtype(:,i_tr);
    
    % which RDK shows event?
    conmat.trials(i_tr).eventRDK = conmat.mats.eventRDK(:,i_tr);
    
    % eventdirection ((according to RDK.event.direction) [1 2 3 4])
    conmat.trials(i_tr).eventdirection = conmat.mats.eventdirection(:,i_tr);
    
    % event onset frames
    conmat.trials(i_tr).event_onset_frames = conmat.mats.event_onset_frames(:,i_tr);
    
    % event onset times
    conmat.trials(i_tr).event_onset_times = conmat.mats.event_onset_times(:,i_tr);
    
    % pre-cue frames
    conmat.trials(i_tr).pre_cue_frames = conmat.mats.pre_cue_frames(:,i_tr);
    
    % pre-cue times
    conmat.trials(i_tr).pre_cue_times = conmat.mats.pre_cue_times(:,i_tr);
    
    % post-cue times
    conmat.trials(i_tr).post_cue_times = p.stim.time_postcue;
    
    % post-cue frames
    conmat.trials(i_tr).post_cue_frames = p.stim.time_postcue*p.scr_refrate;
end



    

end

