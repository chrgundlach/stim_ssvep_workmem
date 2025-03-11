function [timing,key,resp] = pres_SSVEP_workmem(p, ps, key, RDK, conmat, blocknum, flag_training)
% presents experiment SSVEP_FShiftBase
%   p               = parameters
%   ps              = screen parameters
%   RDK             = RDK parameters
%   blocknum        = number of block
%   flag_training   = flag for trainig (1) or experiment (0)


% to do: how to structure trial and trigge
% trialstart -> cue -> post-cue period -> retention interval with no stimuli
% check line 123

%% adaptations for training
if flag_training == 1
    blocknum_present = blocknum;
    blocknum = 1;
end

%% prepare datapixx output
bufferAddress = 8e6;
% % samplesPerTrigger = 1;
triggersPerRefresh = 1; % send one per refresh (i.e. p.scr_refrate at 480 Hz)


%% initialize some RDK settings
% define input for RDK init function
RDKin.scr = ps; RDKin.scr.refrate = p.scr_refrate;
RDKin.Propixx = p.scr_imgmultipl;
RDKin.RDK = RDK;
RDKin.crs = p.crs;

%% create texture of rectangle
% create a white rectangle texture
tx.rectangletex = ones(RDK.RDK(1).dot_size);

% make it a texture
tx.rectTexture = Screen('MakeTexture', ps.window, tx.rectangletex);

%% loop for each trial
trialindex = find(conmat.mats.block==blocknum);

% send start trigger
if flag_training~=1
    PRPX_sendRecTrigger('start')
end

% Wait for release of all keys on keyboard, then sync us to retrace:
KbWait(ps.RespDev,1); % wait for releasing keys on indicated response device

% create keayboard queue
KbQueueCreate(ps.RespDev)

if flag_training~=1
    fprintf(1,'\nexperiment block %1.0f - Praesentation - Trial:', blocknum)
else
    fprintf(1,'\ntraining block %1.0f - Praesentation - Trial:', blocknum_present)
end
ttt=WaitSecs(0.3);

% loop across trials
for i_tr = 1:numel(trialindex)
    fprintf('%1.0f',mod(i_tr,10))
    inittime=GetSecs;
    %% initialize trial structure, RDK, cross, logs
    RDKin.trial = struct('duration',conmat.trials(trialindex(i_tr)).postcue_times+conmat.trials(trialindex(i_tr)).precue_times,...
        'frames',conmat.trials(trialindex(i_tr)).postcue_frames+conmat.trials(trialindex(i_tr)).precue_frames,...
        'cue',conmat.trials(trialindex(i_tr)).precue_frames+1,'event','none');
    RDKin.RDK.event.type = 'none';
    
    % assign RDKs
    t.idx = size(conmat.trials(trialindex(i_tr)).RDK2display);
    RDKin.RDK.RDK = repmat(RDK.RDK(conmat.trials(trialindex(i_tr)).RDK2display(1)),t.idx);
    for i_count = 1:numel(conmat.trials(trialindex(i_tr)).RDK2display)
        RDKin.RDK.RDK(i_count) = RDK.RDK(conmat.trials(trialindex(i_tr)).RDK2display(i_count));
    end
    RDKin.RDK.RDK = RDKin.RDK.RDK';

    RDKin.RDK.presentedRDKs = conmat.trials(trialindex(i_tr)).RDK2display;
    
    % initialize positions and colors and luminance: lummat is the important for texture alpha
    % for stimulus tracking time window (pre-cue and post-cue)
    [colmat,dotmat,dotsize,rdkidx,frames,lummat] = RDK_init_SSVEP_workmem(RDKin.scr,RDKin.Propixx,RDKin.RDK,RDKin.trial,RDKin.crs);

    %%% there is an error in the allocation of the RDK structure!


    % append frames to also include retention interval
    append_frnum = conmat.trials(trialindex(i_tr)).retention_frames/p.scr_imgmultipl;
    frames.flips = frames.flips+append_frnum;
    frames.pertrial = frames.pertrial+conmat.trials(trialindex(i_tr)).retention_frames;
    colmat = cat(3, colmat, repmat(p.scr_color'.*[1 1 1 0]',[1 size(colmat,2) append_frnum]));
    % index last position of all dots (on frame)
    [lastpos_on, lastpos_full] = deal(nan(size(lummat,1),1));
    for i_dot = 1:size(lummat,1)
        lastpos_on(i_dot) = find(lummat(i_dot,:)>0,1,'last'); % find last on position for each dot any nonzero luminance
        lastpos_full(i_dot) = find(lummat(i_dot,:)==1,1,'last'); % find last on position for each dot full on
    end
    append_dotmat = nan([size(dotmat,[1:2]),1]);
    for i_dot = 1:size(dotmat,2)
        append_dotmat(:,i_dot) = dotmat(:,i_dot,lastpos_on(i_dot)); % index position of dot during their last on frame
    end
    dotmat = cat(3,dotmat,repmat(append_dotmat,[1,1,append_frnum]));
    dotsize = cat(2,dotsize,repmat(dotsize(:,end),[1, append_frnum]));
    rdkidx = cat(2,rdkidx,repmat(rdkidx(:,end),[1, append_frnum]));
    lummat = cat(2,lummat,zeros(size(lummat,1), append_frnum));
    
    % append frames to also include first frame of probe stimulus shown (all object on screen)
    frames.flips = frames.flips+1;
    frames.pertrial = frames.pertrial + p.scr_imgmultipl;
    for i_dot = 1:numel(lastpos_on)
        colmat(:,i_dot,frames.flips) = colmat(:,i_dot,lastpos_full(i_dot));
    end
    dotmat(:,:,frames.flips) = dotmat(:,:,frames.flips-1);
    dotsize(:,frames.flips) = dotsize(:,frames.flips-1);
    rdkidx(:,frames.flips) = rdkidx(:,frames.flips-1);
    lummat(:,frames.flips) = ones(size(lummat,1),1);

    % figure; plot(squeeze(dotmat(1,1:5,:))')

    % now get the angles of each of the stimuli
    stimangles = conmat.trials(trialindex(i_tr)).stimangles(unique(conmat.trials(trialindex(i_tr)).RDK2display),:)';
    eventangles = conmat.trials(trialindex(i_tr)).eventangles(unique(conmat.trials(trialindex(i_tr)).RDK2display),:)';
    anglemat = repmat(stimangles(:),p.scr_imgmultipl,frames.flips);
    anglemat(:,end) = repmat(eventangles(:),p.scr_imgmultipl,1);

    % from dotmat to rectmat (required for drawing of rectangles)
    rectmat = nan([4 size(dotmat,2:3)]);
    for i_rect = 1:size(dotmat,2)
        rectmat(:,i_rect,:) = ...
            CenterRectOnPointd([0 0 RDK.RDK(1).dot_size], squeeze(dotmat(1,i_rect,:)),squeeze(dotmat(2,i_rect,:)))';
        % rectmat(:,i_rect,:) = ...
        %     CenterRectOnPointd([1920 1080 RDK.RDK(1).dot_size], squeeze(dotmat(1,i_rect,:)),squeeze(dotmat(2,i_rect,:)))';
    end
    rectmat = rectmat + [ps.xCenter; ps.yCenter; ps.xCenter; ps.yCenter];
    
    % total frames
    totalframes = frames.flips;
       
     % initialize fixation cross 
    crossmat = zeros(1,totalframes); % crossmat has the length of flips, not frames; this means the cue can only change to a full flip and the pre-cue event can only happen to a full flip.
    crossmat(1:conmat.trials(trialindex(i_tr)).precue_frames/p.scr_imgmultipl) = p.FixTex.default;
    if conmat.trials(trialindex(i_tr)).precue_eventnum
        crossmat(conmat.trials(trialindex(i_tr)).precue_event_onset_frames/p.scr_imgmultipl :... % from pre-cue onset frame (divided by 4)
            (conmat.trials(trialindex(i_tr)).precue_event_onset_frames+p.stim.precue_event.length*p.scr_refrate)/p.scr_imgmultipl)... % to end of pre-cue duration (divided by 4)
            = p.FixTex.event(conmat.trials(trialindex(i_tr)).precue_eventid); % write the pre-cue event id
    end
    crossmat(conmat.trials(trialindex(i_tr)).precue_frames/p.scr_imgmultipl+1:end) = p.FixTex.cue(conmat.trials(trialindex(i_tr)).cue);
    
    % preallocate timing
    totaltotalframes = totalframes + max(p.stim.time_response)*p.scr_refrate/p.scr_imgmultipl;
    timing(i_tr) = struct('VBLTimestamp',NaN(1,totaltotalframes),'StimulusOnsetTime',NaN(1,totaltotalframes),...
        'FlipTimestamp',NaN(1,totaltotalframes),'Missed',NaN(1,totaltotalframes));
    
    %% set up responses
    %setup key presses
    key.presses{i_tr}=nan(totaltotalframes,sum(key.keymap));
    key.presses_t{i_tr}=nan(totaltotalframes,sum(key.keymap));
    
    resp(i_tr).trialnumber              = trialindex(i_tr);
    resp(i_tr).blocknumber              = conmat.trials(trialindex(i_tr)).blocknum;
    resp(i_tr).trialtype                = conmat.trials(trialindex(i_tr)).trialtype;
    resp(i_tr).condition                = conmat.trials(trialindex(i_tr)).condition;
    resp(i_tr).RDK2display              = sort(conmat.trials(trialindex(i_tr)).RDK2display(:));
    resp(i_tr).cue                      = conmat.trials(trialindex(i_tr)).cue; % attended left (1) or right (4)
    resp(i_tr).color_attended           = RDK.RDK(resp(i_tr).cue).col(1,:);
    resp(i_tr).freq_attended            = RDK.RDK(resp(i_tr).cue).freq;
    resp(i_tr).cue_onset_fr             = conmat.trials(trialindex(i_tr)).precue_frames + 1;
    resp(i_tr).cue_onset_t_est          = (conmat.trials(trialindex(i_tr)).precue_frames + 1)/p.scr_refrate*1000;
    resp(i_tr).cue_onset_t_meas         = nan; % measured onset time for cue
    resp(i_tr).precue_frames            = conmat.trials(trialindex(i_tr)).precue_frames;
    resp(i_tr).precue_times             = conmat.trials(trialindex(i_tr)).precue_times;
    resp(i_tr).postcue_frames           = conmat.trials(trialindex(i_tr)).postcue_frames;
    resp(i_tr).postcue_times            = conmat.trials(trialindex(i_tr)).postcue_times;
    resp(i_tr).retention_frames         = conmat.trials(trialindex(i_tr)).retention_frames;
    resp(i_tr).retention_times          = conmat.trials(trialindex(i_tr)).retention_times;
    resp(i_tr).eventtype                = conmat.trials(trialindex(i_tr)).eventtype; % 1 = match; 2 = non-match
    resp(i_tr).eventRDK                 = conmat.trials(trialindex(i_tr)).eventRDK;
    resp(i_tr).eventcolor               = RDK.RDK(resp(i_tr).eventRDK).col(1,:);
    resp(i_tr).eventfreq                = RDK.RDK(resp(i_tr).eventRDK).freq(1,:);
    resp(i_tr).stimangles               = conmat.trials(trialindex(i_tr)).stimangles(resp(i_tr).RDK2display,:);
    resp(i_tr).eventangles              = conmat.trials(trialindex(i_tr)).eventangles(resp(i_tr).RDK2display,:);
    resp(i_tr).precue_eventnum          = conmat.trials(trialindex(i_tr)).precue_eventnum;
    resp(i_tr).precue_eventid           = conmat.trials(trialindex(i_tr)).precue_eventid;
    resp(i_tr).precue_eventtype         = conmat.trials(trialindex(i_tr)).precue_eventtype;
    resp(i_tr).precue_event_onset_fr    = conmat.trials(trialindex(i_tr)).precue_event_onset_frames;
    resp(i_tr).precue_event_onset_t_est = conmat.trials(trialindex(i_tr)).precue_event_onset_times;
    
    %% set up datapixx trigger vector
    % prepare datapixx scheduler
    doutWave = zeros(1,length(crossmat)*p.scr_imgmultipl);
    
    % write trigger numbers into doutwave
    % trial start trigger
    doutWave(1) = p.trig.tr_start;

    % add condition trigger when cue is presented
    % condition trigger for cue
    resp(i_tr).triggernum_cue = ...
        p.trig.tr_con_cue(resp(i_tr).condition); % condition 1 2 3 4 5 6
    doutWave(resp(i_tr).cue_onset_fr) = resp(i_tr).triggernum_cue;
    
    % condition trigger for retention
    resp(i_tr).triggernum_retent = ...
        p.trig.tr_con_retent(resp(i_tr).condition); % condition 1 2 3 4 5 6
    doutWave(resp(i_tr).cue_onset_fr + resp(i_tr).postcue_frames) = resp(i_tr).triggernum_retent;
    
    % condition trigger for event
    resp(i_tr).triggernum_event = ...
        p.trig.tr_con_event(resp(i_tr).eventtype, resp(i_tr).condition); % condition 1 2 3 4 5 6
    doutWave(resp(i_tr).cue_onset_fr + resp(i_tr).postcue_frames + resp(i_tr).retention_frames) = resp(i_tr).triggernum_event;
    
    % add pre-event trigger
    if resp(i_tr).precue_eventnum > 0
        resp(i_tr).triggernum_precue_event = ...
            p.trig.tr_precue_event(resp(i_tr).precue_eventtype);
        doutWave(resp(i_tr).precue_event_onset_fr) = resp(i_tr).triggernum_precue_event;
    end
     
    doutWave = [doutWave;zeros(triggersPerRefresh-1,numel(doutWave))]; doutWave=doutWave(:);
    samplesPerFlip = triggersPerRefresh * p.scr_imgmultipl;
    % figure; plot(doutWave)       
    
    %% getting wait ITI time and then start drawing routines
    % wait
    if i_tr > 1
        crttime2 = GetSecs;
        t.timetowait = (p.ITI(1)/1000)-(crttime2 - crttime);
        ttt=WaitSecs(t.timetowait);
    end

    %%%%
    % include waiting routine here!
    %%%%
    % draw fixation cross
    Screen('DrawTextures', ps.window, p.FixTex.default, [], p.crs.rects);
    
    % send 0 before again to reset everything
    Datapixx('SetDoutValues', 0);
    Datapixx('RegWrRd');
    
    % write outsignal
    Datapixx('WriteDoutBuffer', doutWave, bufferAddress);
    % disp(Datapixx('GetDoutStatus'));
    Datapixx('SetDoutSchedule', 0, [samplesPerFlip, 2], numel(doutWave), bufferAddress); % 0 - scheduleOnset delay, [samplesPerFlip, 2] - sheduleRate in samples/video frame, framesPerTrial - maxSheduleFrames
    Datapixx('StartDoutSchedule');

    %% keyboard
    % start listening to keyboard
    KbQueueStart(ps.RespDev);
    KbQueueFlush(ps.RespDev);
    
    % flip to get everything synced
    Screen('Flip', ps.window, 0);
    
    %% loop across frames for pre-cue | tracking | retention | first frame of probe
    for i_fl = 1:frames.flips
        %% Drawing
        % draw rectangles
        Screen('DrawTextures', ...  % draw rectangles
            ps.window, ...          % in specified window
            tx.rectTexture, ...     % which texture
            [], ...                 % source rectangle (not specified)
            rectmat(:,:,i_fl),...   % rectangle position
            anglemat(:,i_fl)', ...  % rotation angles
            0, ...                  % filter mode for texture display (if resized)
            [], ...                 % global alpha not specified
            colmat(:,:,i_fl) ...    % color matrix
            );  

        % fixation cross
        Screen('DrawTextures', ps.window, crossmat(i_fl), [], p.crs.rects);
    %     Screen('Flip', ps.window, 0);
    % end
        
        %% start trigger schedule and start listening to response device
        if i_fl == 1 % send the trigger with the start of the 1st flip
            Datapixx('RegWrVideoSync');
        end
        
        % Flip
        [timing(i_tr).VBLTimestamp(i_fl), timing(i_tr).StimulusOnsetTime(i_fl), timing(i_tr).FlipTimestamp(i_fl), timing(i_tr).Missed(i_fl)] = Screen('Flip', ps.window, 0);

        % get image
        %current_display = Screen('GetImage',ps.window);
        %imwrite(current_display,'./stims/all_nocue_centralredblue_periredgreen.png')
        
        % send trigger/save timing/ reset timing
        if i_fl == 1
            % start trigger
            starttime=GetSecs;
            KbEventFlush(ps.RespDev); % flush keyboard
        elseif i_fl == frames.flips
            starttime_probe = GetSecs; % onset time for probe display
        end
        
        %% check for button presses
        [key.pressed, key.firstPress]=KbQueueCheck(ps.RespDev);
        key.presses{i_tr}(i_fl,:)=key.firstPress(key.keymap)>1;
        key.presses_t{i_tr}(i_fl,:)=(key.firstPress(key.keymap)-starttime).*key.presses{i_tr}(i_fl,:);
        if any(key.firstPress(key.keymap)>1)
            lptwrite(1,find(key.firstPress(key.keymap),1,'first'),500);
        end

    end
    
    %% probe presentation until response or maximum time window
    i_fl_r = i_fl+1;
    crttime3 = GetSecs;
    while ~any(key.firstPress(key.keymap)>1) & crttime3 - starttime_probe <= max(p.stim.time_response)
        %% display probe array
        % draw rectangles
        Screen('DrawTextures', ...  % draw rectangles
            ps.window, ...          % in specified window
            tx.rectTexture, ...     % which texture
            [], ...                 % source rectangle (not specified)
            rectmat(:,:,end),...   % rectangle position
            anglemat(:,end)', ...  % rotation angles
            0, ...                  % filter mode for texture display (if resized)
            [], ...                 % global alpha not specified
            colmat(:,:,end) ...    % color matrix
            );  

        % fixation cross
        Screen('DrawTextures', ps.window, crossmat(end), [], p.crs.rects);

        % Flip
        [timing(i_tr).VBLTimestamp(i_fl_r), timing(i_tr).StimulusOnsetTime(i_fl_r), timing(i_tr).FlipTimestamp(i_fl_r), timing(i_tr).Missed(i_fl_r)] = ...
            Screen('Flip', ps.window, 0);


        %% check for button presses
        [key.pressed, key.firstPress]=KbQueueCheck(ps.RespDev);
        key.presses{i_tr}(i_fl_r,:)=key.firstPress(key.keymap)>1;
        key.presses_t{i_tr}(i_fl_r,:)=(key.firstPress(key.keymap)-starttime).*key.presses{i_tr}(i_fl_r,:);
        if any(key.firstPress(key.keymap)>1)
            lptwrite(1,find(key.firstPress(key.keymap),1,'first'),500);
        end
        
        
        % update running flip
        i_fl_r = i_fl_r+1;
        crttime3 = GetSecs;
    end

    
    
    %% ITI
    % fixation cross
    Screen('DrawTextures', ps.window, crossmat(end), [], p.crs.rects);
    
    % flip
    Screen('Flip', ps.window, 0);
    
    % get time for end of trial start of ITI
    crttime = GetSecs;
    
    %%%%%%%%
    % do behavioral calculation
    key.presses{i_tr}(1,:)=[];
    key.presses_t{i_tr}(1,:)=[];

    key.presses{i_tr}(find(isnan(key.presses{i_tr}),1,'first'):end,:)=[];
    key.presses_t{i_tr}(find(isnan(key.presses_t{i_tr}),1,'first'):end,:)=[];
    
    % get frame and timing of button press onsets
    resp(i_tr).button_presses_fr=nan(max(sum(key.presses{i_tr})),size(key.presses{i_tr},2));
    resp(i_tr).button_presses_t=resp(i_tr).button_presses_fr;
    for i_bt = 1:size(key.presses{i_tr},2)
        try
            resp(i_tr).button_presses_fr(1:sum(key.presses{i_tr}(:,i_bt)),i_bt)=...
                find(key.presses{i_tr}(:,i_bt));
            resp(i_tr).button_presses_t(1:sum(key.presses{i_tr}(:,i_bt)),i_bt)=...
                key.presses_t{i_tr}(find(key.presses{i_tr}(:,i_bt)),i_bt)*1000; % in ms
        catch
            resp(i_tr).button_presses_fr(:,i_bt)=nan;
        end
    end
    
    % check for hits {'hit','miss','error'}
    resp(i_tr).event_response_type = {}; %{'hit','miss','error'}
    resp(i_tr).event_response_RT = []; % reaction time in ms (after event or closest to other event)
    % all relevant presses 
    t.presses = resp(i_tr).button_presses_t(:,key.keymap_ind==key.SAME | key.keymap_ind==key.diff);
    % first define response windows
    if any(~isnan(resp(i_tr).eventtype))
        t.respwin = ((resp(i_tr).precue_times+resp(i_tr).postcue_times+resp(i_tr).retention_times)+(p.targ_respwin/1000))*1000;
    end
    
    % check for button presses in time window
    t.idx = t.presses>=t.respwin(:,1) & t.presses<=t.respwin(:,2);
    % not in any response window? --> miss; no reaction time
    if ~any(t.idx,"all")
        resp(i_tr).event_response_type = 'miss';
        % no reaction time
        resp(i_tr).event_response_RT = nan;
    else
        if find(any(t.idx,1)) == resp(i_tr).eventtype % check for button presses of same for same trials and diff for diff trials
            resp(i_tr).event_response_type = 'hit';
        else
            resp(i_tr).event_response_type = 'error';
        end

        % reaction time?
        resp(i_tr).event_response_RT = t.presses(t.idx)- ...
            (resp(i_tr).precue_times + resp(i_tr).postcue_times + resp(i_tr).retention_times)*1000;
    end

    % check for pre-cue events
    % check for hits {'hit','miss','error'}
    resp(i_tr).precue_event_response_type = {}; %{'hit','miss','error'}
    resp(i_tr).precue_event_response_RT = []; % reaction time in ms (after event or closest to other event)

    % first define response windows
    if any(~isnan(resp(i_tr).precue_eventtype))
        t.respwin = (resp(i_tr).precue_event_onset_t_est+(p.targ_respwin/1000))*1000;
    end

    % check for button presses in time window
    t.idx = find(t.presses>=t.respwin(:,1) & t.presses<=t.respwin(:,2));
    
    % not in any response window? --> miss; no reaction time
    if resp(i_tr).precue_eventnum < 1
        resp(i_tr).precue_event_response_type = 'noevent';
        % no reaction time
        resp(i_tr).precue_event_response_RT = nan;
    elseif ~any(t.idx,"all") & resp(i_tr).precue_eventnum == 1
        resp(i_tr).precue_event_response_type = 'miss';
        % no reaction time
        resp(i_tr).precue_event_response_RT = nan;
    else
        % check for first respons (in case of two presses: should be rare but might happen)
        [~, t.idx2] = min(t.presses(t.idx));
        t.idx3 = false(size(t.presses)); t.idx3(t.idx(t.idx2))=true;

        if find(any(t.idx3,1)) == resp(i_tr).precue_eventtype % check for button presses of same for same trials and diff for diff trials
            resp(i_tr).precue_event_response_type = 'hit';
        else
            resp(i_tr).precue_event_response_type = 'error';
        end

        % reaction time?
        resp(i_tr).precue_event_response_RT = t.presses(t.idx3)- ...
            resp(i_tr).precue_event_onset_t_est*1000;
    end
    
end

% stop response recording
KbQueueRelease(ps.RespDev);

% send start trigger
if flag_training~=1
    PRPX_sendRecTrigger('stop')
end

% close textures
Screen('Close', tx.rectTexture);



end