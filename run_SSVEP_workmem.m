function [] = run_SSVEP_workmem(sub,flag_training, flag_block)
% run_SSVEP_workmem(sub,flag_training, flag_block)
%   runs experiment SSVEP_FShiftBase
%       sub:            participant number
%       flag_training:  1 = do training
%       flag_block:     1 = start with block 1
%           e.g. run_SSVEP_workmemr(1,1,1)
% 
%   study looking into working memory and SSVEP relationship
%
% to be decided?
%   - stimulus sizes
%       all for three stimuli
%       [ext_width ext_height]=visualangle([50 15],[120],[1920 1080]./2,[63 36])  1.56° x 0.48° --> favorite size
%       [ext_width ext_height]=visualangle([64 16],[120],[1920 1080]./2,[63 36])  2° x 0.5°
%       [ext_width ext_height]=visualangle([308 308],[120],[1920 1080]./2,[63 36])  9.5° x 9.7°
%       [ext_width ext_height]=visualangle([290 290],[120],[1920 1080]./2,[63 36])  9.01° x 9.15°
%       [ext_width ext_height]=visualangle([260 290],[120],[1920 1080]./2,[63 36])  8.09° x 9.15°--> favorite area
%       [ext_width ext_height]=visualangle([-192],[120],[1920 1080]./2,[63 36])  -6.0°
%       [ext_width ext_height]=visualangle([-160],[120],[1920 1080]./2,[63 36])  -5.0°          --> favorite shift
%
%       all for three stimuli
%       [ext_width ext_height]=visualangle([50 15],[120],[1920 1080]./2,[63 36])  1.56° x 0.48° --> favorite size
%       [ext_width ext_height]=visualangle([64 16],[120],[1920 1080]./2,[63 36])  2° x 0.5°
%       [ext_width ext_height]=visualangle([160 238],[120],[1920 1080]./2,[63 36])  5° x 7.5°
%       [ext_width ext_height]=visualangle([176 254],[120],[1920 1080]./2,[63 36])  5.5° x 8°
%       [ext_width ext_height]=visualangle([176 270],[120],[1920 1080]./2,[63 36])  5.5° x 8.5°
%       [ext_width ext_height]=visualangle([192 270],[120],[1920 1080]./2,[63 36])  6° x 8.5°   --> favorite shift
%       [ext_width ext_height]=visualangle([-120],[120],[1920 1080]./2,[63 36])  -3.5°          
%       [ext_width ext_height]=visualangle([-128],[120],[1920 1080]./2,[63 36])  -3.5°          --> favorite shift
%
%       original study Vogel et al.
%           RDK region rectangular 4° x 7.3° shifted left and right by 3°, stimuli 0.65° (seems small); two stimuli or four stimuli

% Christopher Gundlach,  Leipzig, 2025

if nargin < 3
    help run_SSVEP_workmem
    return
end

%% parameters
% sub = 1; flag_training = 1; flag_block = 1;
% design
p.sub                   = sub;                  % subject number
p.flag_block            = flag_block;           % block number to start
p.flag_training         = flag_training;        % do training

p.ITI                   = [1000 1000];          % inter trial interval in ms
p.targ_respwin_main     = [200 1500];           % time window for responses in ms
p.targ_respwin_precue   = [200 1000];           % time window for responses in ms


% screen
p.scr_num               = 1;                    % screen number
p.scr_res               = [1920 1080];          % resolution
p.scr_refrate           = 480;                  % refresh rate in Hz (e.g. 85)
p.scr_color             = [0.05 0.05 0.05 1];      % default: [0.05 0.05 0.05 1]; ; color of screen [R G B Alpha]
p.scr_imgmultipl        = 4;

p.isol.override         = [0.4706 0.1882 0 1; 0 0.3498 0.8745 1;0 0.4392 0 1]; % these are the ones used for p.isol.bckgr            = p.scr_color(1:3)+0.2;



% stimplan
p.stim.condition        = [1 2 3 4 5 6];    
                        % [cue_left_3t3d; cue_right_3t3d; cue_left_6t; cue_right_6t; cue_left_3t; cue_right_3t]
                        % [RDK1 RDK3; RDK4 RDK6] [RDK1 RDK3; RDK4 RDK6] [RDK1 RDK2; RDK4 RDK5] [RDK1 RDK2; RDK4 RDK5] [RDK1; RDK4] [RDK1; RDK4] 
p.stim.con_targstim     = [3 3 6 6 3 3];
p.stim.con_diststim     = [3 3 0 0 0 0];
p.stim.con_repeats      = [85 85 85 85 85 85];  % trial number/repeats for each eventnum and condition
p.stim.con_repeats_catch= [15 15 15 15 15 15];  % trial number/repeats for each eventnum and condition
p.stim.trialnum_train   = [12];
p.stim.time_postcuetrack= [1.5 2];          % post.cue tracking time in s [upper lower]
p.stim.time_postcuetrack_catch = [0.7 1.5]; % post.cue tracking time in s [upper lower] for catch trials
p.stim.time_precue      = [1.5 2];          % precue time in s; [upper lower] for randomization
p.stim.time_retention   = [1 1];            % time for retention period time in s; [upper lower] for randomization
p.stim.time_response    = [1 2];            % time range of response time window [min max]
p.stim.event.ratio      = [1 1];            % ratio of events; probed display matches vs doesn't match
p.stim.angles           = [0 45 90 135];    % possible rotations of squares
p.stim.event.rotations  = [-45 45];         % possible rotations of squares for events
p.stim.blocknum         = 20;               % number of blocks
p.stim.ITI              = [1 1];            % ITI range in seconds

% pre-cue events
p.stim.precue_event.num         = [0 0 0 1]; % ratio of no precue-events (0) and precue-events(1); set 0 to turn precue-events off
p.stim.precue_event.targets     = [1 3];     % defines which are the target events: 1 - left shorter, 2 - left longer, 3 - right shorter, 4 - right longer
p.stim.precue_event.length      = 0.15;      % length of precue-event in s              
p.stim.precue_event.min_onset   = 0.3;       % min time before precue-event onset in s
p.stim.precue_event.min_offset  = 0.4;       % min offset from precue-event end to end of trial in s


% introduce RDK structure
% RDK.RDK(1).size         = [290 290];                    % width and height of RDK in pixel; only even values [308 = 9.6°] [290 = 9°]
% RDK.RDK(1).size         = [260 290];                    % width and height of RDK in pixel; only even values [8.1° x 9.1°]
RDK.RDK(1).size         = [192 270];                    % width and height of RDK in pixel; only even values [5° x 7.5°]
RDK.RDK(1).centershift  = [0 0];                        % position of RDK center; x and y deviation from center in pixel
RDK.RDK(1).col          = [1 1 1 1; p.scr_color(1:3) 0];% "on" and "off" color [assigned later]
RDK.RDK(1).freq         = 0;                            % flicker frequency, frequency of a full "on"-"off"-cycle
RDK.RDK(1).mov_freq     = 120;                          % Defines how frequently the dot position is updated; 0 will adjust the update-frequency to your flicker frequency (i.e. dot position will be updated with every "on"-and every "off"-frame); 120 will update the position for every frame for 120Hz or for every 1. quadrant for 480Hz 
% RDK.RDK(1).num          = 3;                            % number of dots % 85
RDK.RDK(1).num          = 2;                            % number of dots % 85
RDK.RDK(1).mov_speed    = 2;                            % movement speed in pixel
RDK.RDK(1).mov_dir      = [0 1; 0 -1; -1 0; 1 0];       % movement direction  [0 1; 0 -1; -1 0; 1 0] = up, down, left, right
RDK.RDK(1).dot_size     = [50 15];                      % size of rectangles
RDK.RDK(1).shape        = 1;                            % 1 = square RDK; 0 = ellipse/circle RDK;
RDK.RDK(1).id           = 1;                            % identifier for RDK


% p.stim.pos_shift        = [-192 0; 192 0];              % position shift in pixel for stimuli in periphery [255 = 7.8°] either left or right
p.stim.pos_shift        = [-128 0; 128 0];              % position shift in pixel for stimuli in periphery [255 = 7.8°] either left or right

p.stim.freqs            = [17 20 23 26];                % possible frequencies of different RDKs
p.stim.colors           = ...                           % "on" and "off" color
    {[0.4706 0.1882 0 1; p.scr_color(1:3) 0];...
    [0 0.3498 0.8745 1; p.scr_color(1:3) 0]};
    % plot_colorwheel([0.4706 0.1882 0; 0 0.3498 0.8745],'ColorSpace','propixxrgb','LAB_L',50,'NumSegments',60,'AlphaColWheel',1,'LumBackground',100, 'disp_LAB_vals', 1)
p.stim.color_names      = {'redish';'blue'};
p.stim.RDKsamecolors    = [1 1 2 1 1 2]; % which ones of the RDKs share the same color?
p.stim.RDKsamefreqs     = [1 2 2 3 4 4];
p.stim.RDKpos_shifts    = [1 1 1 2 2 2];
 
% fixation cross
p.crs.color             = [0.4 0.4 0.4 1];      % color of fixation cross
p.crs.size              = 12;                   % size of fixation
p.crs.width             = 2;                    % width of fixation cross
p.crs.cutout            = 0;                    % 1 = no dots close to fixation cross
p.crs.event             = 1;                    % pixels that disapear for event

% trigger
p.trig.rec_start        = 253;                  % trigger to start recording
p.trig.rec_stop         = 254;                  % trigger to stop recording
p.trig.tr_start         = 77;                   % trial start; main experiment
p.trig.tr_con_cue       = [1 2 3 4 5 6 ]*10;    % condition type for cue presentation
p.trig.tr_con_retent    = p.trig.tr_con_cue+1;  % condition type for start of retention period
p.trig.tr_con_event     = p.trig.tr_con_cue+5;  % condition type for start of retention period | match
p.trig.tr_con_event(2,:)= p.trig.tr_con_cue+6;  % condition type for start of retention period | change
p.trig.tr_precue_event  = [81 82];              % pre-cue event target or distractor
p.trig.button           = [155 156];            % button press (match vs nonmatch)

% possible condition triggers:
% {[1 101 201 111 121 211 221]; [2 102 202 112 122 212 222]; [3 103 203 113 123 213 223]; ...
% [4 104 204 114 124 214 224]; [5 105 205 115 125 215 225]; [6 106 206 116 126 216 226]}

% logfiles
p.log.path              = '/home/stimulation120/matlab/user/christopher/stim_ssvep_workmem/logfiles/';
p.log.exp_name          = 'SSVEP_workmem';
p.log.add               = '_a';


%% check for logfile being present
filecheck=dir(sprintf('%sVP%02.0f_timing*',p.log.path,p.sub));
if ~isempty(filecheck)
    reply = input(sprintf('\nVP%02.0f existiert bereits. Datei überschreiben? [j/n]... ',p.sub),'s');
    if strcmp(reply,'j')
        p.filename = sprintf('VP%02.0f_timing',p.sub);
    else
        [~, name_ind]=max(cellfun(@(x) numel(x), {filecheck.name}));
        p.filename = sprintf('%s%s',filecheck(name_ind).name(1:end-4),p.log.add);
    end
else
    p.filename = sprintf('VP%02.0f_timing',p.sub);
end


%% Screen init
ps.input = struct('ScrNum',p.scr_num,'RefRate',p.scr_refrate,'PRPXres',p.scr_res,'BckGrCol',p.scr_color,'PRPXmode',2);
[~, ps.screensize, ps.xCenter, ps.yCenter, ps.window, ps.framerate, ps.RespDev, ps.keymap] = PTExpInit_GLSL(ps.input,1);
% ps.xCenter = 1920/2; ps.yCenter=1080/2;

% some initial calculations
% fixation cross
ps.center = [ps.xCenter ps.yCenter];
p.crs.half = p.crs.size/2;
p.crs.bars = [-p.crs.half p.crs.half 0 0; 0 0 -p.crs.half p.crs.half];

% shift into 4 quadrants (running with 480 Hz)
ps.shift = [-ps.xCenter/2, -ps.yCenter/2; ps.xCenter/2, -ps.yCenter/2;... % shifts to four quadrants: upper left, upper right, lower left, lower right
    -ps.xCenter/2, ps.yCenter/2; ps.xCenter/2, ps.yCenter/2];

p.crs.rects = zeros(4,p.scr_imgmultipl);
baserect = [0 0 p.crs.size p.crs.size*2 + p.crs.width-1];
if p.scr_imgmultipl > 1
    for rect = 1:p.scr_imgmultipl
        p.crs.rects(:,rect) = CenterRectOnPoint(baserect,ps.xCenter+ps.shift(rect,1),ps.yCenter+ps.shift(rect,2));
    end
else
    p.crs.rects = CenterRectOnPoint(baserect, ps.xCenter, ps.yCenter);
end

%% keyboard and ports setup ???
KbName('UnifyKeyNames')
Buttons = [KbName('ESCAPE') KbName('Q') KbName('SPACE') KbName('j') KbName('n') KbName('s') KbName('d')];
RestrictKeysForKbCheck(Buttons);
key.keymap=false(1,256);
key.keymap(Buttons) = true;
key.keymap_ind = find(key.keymap);
[key.ESC, key.SECRET, key.SPACE, key.YES, key.NO, key.SAME, key.diff] = deal(...
    Buttons(1),Buttons(2),Buttons(3),Buttons(4),Buttons(5),Buttons(6),Buttons(7));

%% start experiment
% initialize randomization of stimulation frequencies and RDK colors
% inititalize RDKs [RDK1 and RDK2 task relevant at center;  RDK3 RDK4 RDK5 not and in periphery]

% rand('state',p.sub)
rng(p.sub,'v4')
% pre-allocate new RDKs
RDK.RDK(2:numel(p.stim.RDKsamecolors)) = deal(RDK.RDK(1));

% randomize colors
t.rndselidx = randsample(numel(p.stim.colors),numel(p.stim.colors));
for i_col = 1:numel(p.stim.colors)
    [RDK.RDK(p.stim.RDKsamecolors==i_col).col] = deal(p.stim.colors{t.rndselidx(i_col)});
end

% randomize frequencies
t.rndselidx = randsample(numel(p.stim.freqs),numel(p.stim.freqs));
for i_freq = 1:numel(p.stim.freqs)
    [RDK.RDK(p.stim.RDKsamefreqs==i_freq).freq] = deal(p.stim.freqs(t.rndselidx(i_freq)));
end

% assign positions
for i_pos = 1:size(p.stim.pos_shift,1)
    [RDK.RDK(p.stim.RDKpos_shifts==i_pos).centershift] = deal(p.stim.pos_shift(i_pos,:));
end

% assign id
for i_rdk = 1:numel(RDK.RDK)
    RDK.RDK(i_rdk).id = i_rdk;
end

% initialize blank variables
timing = []; button_presses = []; resp = []; randmat = [];

% make fixation "cross" textures
[p.FixTex] = MakeFixationTextures(p,ps,RDK);

%% initial training
if p.flag_training
    fprintf(1,'\nTraing starten mit q')
    inp.prompt_check = 0;
    while inp.prompt_check == 0             % loop to check for correct input
        [key.keyisdown,key.secs,key.keycode] = KbCheck;
        if key.keycode(key.SECRET)==1
            flag_trainend = 0; inp.prompt_check = 1;
        end
        Screen('Flip', ps.window, 0);
    end
    
    
    i_bl = 1;
    flag_trainend = 0;
    while flag_trainend == 0 % do training until ended
        %rand('state',p.sub*i_bl) % determine randstate
        rng(p.sub*i_bl,'v4')
        randmat.training{i_bl} = rand_SSVEP_workmem(p, RDK,  1);
        [timing.training{i_bl},button_presses.training{i_bl},resp.training{i_bl}] = ...
            pres_SSVEP_workmem(p, ps, key, RDK, randmat.training{i_bl}, i_bl,1);
        save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
        pres_feedback(resp.training{i_bl},p,ps, key,RDK)
               
        % loop for training to be repeated
        fprintf(1,'\nTraing wiederholen? (j/n)')
        inp.prompt_check = 0;
        while inp.prompt_check == 0             % loop to check for correct input
            [key.keyisdown,key.secs,key.keycode] = KbCheck; 
            if key.keycode(key.YES)==1
                i_bl = i_bl + 1; flag_trainend = 0; inp.prompt_check = 1;
            elseif key.keycode(key.NO)==1
                flag_trainend = 1; inp.prompt_check = 1;
            end
            Screen('Flip', ps.window, 0);
        end  
        
    end
end


%% present each block
% randomization
% rand('state',p.sub);                         % determine randstate
rng(p.sub,'v4')
randmat.experiment = rand_SSVEP_workmem(p, RDK,  0);    % randomization
for i_bl = p.flag_block:p.stim.blocknum
    % start experiment
    [timing.experiment{i_bl},button_presses.experiment{i_bl},resp.experiment{i_bl}] = ...
        pres_SSVEP_workmem(p, ps, key, RDK, randmat.experiment, i_bl,0);
    % save logfiles
    save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
          
    pres_feedback(resp.experiment{i_bl},p,ps, key, RDK)    
end

fprintf(1,'\n\nENDE\n')

%Close everything
Datapixx('SetPropixxDlpSequenceProgram', 0);
Datapixx('RegWrRd');
Datapixx('close');
ppdev_mex('Close', 1);
ListenChar(0);
sca;


end

