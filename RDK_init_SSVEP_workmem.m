function [colmat,dotmat,dotsize,rdkidx,frames,varargout] = RDK_init_SSVEP_workmem(scr,Propixx,RDK,trial,cross)
%RDK_init initializes parameters for RDKs
% This function initializes dot colors and dot positions for every frame.
% The different RDKs are treated as "one big dot cloud". Dots are
% randomized, to ensure that dots of lastly drawn RDK are not systematically
% overlapping the others.
% requires
%       scr.        refrate
%                   shift           shift to 4 quadrants (specific to Propixx)
%       Propixx     4               specific fro Propixx to draw in QUAD4X mode
%       RDK.RDK(x).                 input structure ith RDKs
%                   size            [x y] of ellipse in pixel
%                   centershift     displacement from center in pixel
%                   col             [R G B A; R G B A]; on off colors
%                   mov_freq        frequency for which movement is updated
%                   mov_dir         movement direction [x y; x y;...] in pixel
%                   mov_speed       movement speed in pixel
%                   num             number of dots for each RDK
%                   dot_size        size of dots
%                   flicker_type    'SSVEP'; 'BRBF' SSVEP vs broad band flicker; default: SSVEP
%       RDK.event.  class           'globalmotion'; 'colorchange'
%                   duration        duration in s
%                   coherence       amount of dots showing coherent motion
%                   direction       direction of event movements [x y; x y;...] in pixel for class 'globalmotion'
%                   color           [R G B A; R G B A]; on off color of events for class 'colorchange'
%       trial.      duration        in s
%                   frames
%                   cue             cue frame (time point for syncing)
%                   event.onset     onset frame of event
%                   event.direction direction of event movements [1 2 3 4] in pixel
%                   event.RDK       RDK for event
%       cross.      cutout          no movements close to cross? 1 = no movements
%                   size            size of central fixation cross
%
%       example:
%                   - load('RDK_Init_example_input.mat')
%                   - [colmat,dotmat,dotsize,rdkidx,frames] = RDK_init(RDKin.scr,RDKin.Propixx,RDKin.RDK,RDKin.trial,RDKin.crs);
%                   - [colmat,dotmat,dotsize,rdkidx,frames] = RDK_init(RDKin_adj.scr,RDKin_adj.Propixx,RDKin_adj.RDK,RDKin_adj.trial,RDKin_adj.crs);

% M. Dotzer, C. Gundlach [2017, 2018, 2023, 2025]
% v20231012:
%       - don't allow stimuli to be presented overlappingly
%       - differs from the original function
% v20231012:
%       - includes changes of color as targets
%
%
% To do:

% check for frames
if ~isfield(trial, 'frames')
    frames.pertrial = trial.duration*scr.refrate;
else
    frames.pertrial = trial.frames;
end
frames.flips = frames.pertrial/Propixx; % There are 4 frames per flip in the QUAD4X mode

% check for flicker type and setup default (SSVEP)
if ~isfield(RDK.RDK, 'flicker_type')
    for i_rdk = 1:numel(RDK.RDK)
        RDK.RDK(i_rdk).flicker_type = 'SSVEP';
    end
end

colmat  = [];
dotmat  = [];
dotsize = [];
rdkidx  = [];
lummat  = [];

%% Do computations for every RDK

%% for color and flicker information
for r = 1:numel(RDK.RDK)

    % do classical SSVEP or broadband flicker stimulation?
    switch RDK.RDK(r).flicker_type
        case 'SSVEP'

            % define on frames
            onframes = scr.refrate/RDK.RDK(r).freq/2; %number of "on"-frames

            % setup color matrix sync on cue
            if ~isfield(trial, 'cue') % check for at certain frame (shifting design) vs cued befor trial design
                fr_cue = 1;
            else
                if trial.cue < 2; fr_cue = 1;
                else
                    fr_cue = trial.cue;
                end
            end

            % find number of decimals for frequency to create vector that ends with full cycle (for frequencies with decimal point)
            if mod(RDK.RDK(r).freq,1) > 0
                t.temp = strsplit(num2str(RDK.RDK(r).freq),'.');
                t.dec = length(t.temp{2});
            else
                t.dec = 0;
            end
            % create vector of x*2*10^t.dec seconds that contains two parts ending with a full cycle and that it is at least twice as long
            % as the pre-cue or post-cue window
            if any(scr.refrate*10^t.dec < [frames.pertrial-fr_cue+1 fr_cue-1])
                t.framenum = 2 * scr.refrate * 10^t.dec * ...
                    max(ceil([frames.pertrial-fr_cue+1 fr_cue-1]/scr.refrate*10^t.dec));
            else
                t.framenum = 2*scr.refrate*10^t.dec;
            end
            % create frame vector for MakeFlicker function
            t.frames = 0:t.framenum-1;
            % create raw weights with interpolated flanks (is longer than it should be
            weight_raw = MakeFlicker(t.frames,scr.refrate/RDK.RDK(r).freq);
            % figure; plot(t.frames./scr.refrate, weight_raw)
            % prune weights according to number of pre-cue and post-cue frames
            t.idx = [ (1:fr_cue-1) + numel(weight_raw)/2 - (fr_cue-1)... % pre-cue
                (numel(weight_raw)/2)+(1:frames.pertrial-fr_cue+1)]; % post-cue
            weight = weight_raw(t.idx);
            % figure; plot(t.frames(t.idx)./scr.refrate, weight)
            % figure; plot( weight)

        case 'BRBF'
            % create phase modulated signal
            brbf.amprange = [0 1]; % signal range
            brbf.phase_shift_size =    1; % default 1; 2 increases the amount of phase shift (and spread of frequencies)
            brbf.phase_shift_freqnum = 3; % default = 3; increases the number of added frequencies for phase shift (adds spred of frequecies)
            brbf.time = (1:trial.frames)./scr.refrate;

            weight = (sin( ...
                2 * pi * brbf.time * mean(RDK.RDK(r).freq) + ...
                sum(cell2mat(arrayfun(@(i) sin(2 * pi * (brbf.time + rand(1)) * (i +rand(1)/4 + 1)) .* brbf.phase_shift_size, 1:0.5:brbf.phase_shift_freqnum, 'UniformOutput', false)'),1)...
                ) + (1 + brbf.amprange(1)))./(2/diff(brbf.amprange));

            % test plotting
            % figure; plot(weight)

            %             t.fft_amp = squeeze(abs(fft(detrend(weight),10000,2))*2/size(weight,1));
            %             t.fft_xscale = ((0:size(t.fft_amp,2)-1)/size(t.fft_amp,2)) * scr.refrate;
            %             plot(t.fft_xscale,t.fft_amp(:,:))
            %             title('FFT derived frequency spectra of SSVEP signal')
            %             xlabel('frequency in Hz')
            %             ylabel('amplitude (a.u.)')
            %             xlim([0 100])

    end

    % adjust color values by weight
    % with color targets created if subclass 'colorchange' is defined
    colmat_r = zeros(size(RDK.RDK(r).col,2), RDK.RDK(r).num, frames.pertrial); %setup color matrix: RGB x dot number x frames
    % index which dots are moved coherently
    if ~strcmp(RDK.event.type, 'none')
        t.cohidx = randsample(RDK.RDK(r).num,ceil(RDK.RDK(r).num*RDK.event.coherence));
        t.cohidx2 = false(1,RDK.RDK(r).num); t.cohidx2(t.cohidx)=true;
    end
    for w=1:length(weight)
        if  strcmp(RDK.event.type, 'colorchange') & ... % is color change the event class?
                (w >= trial.event.onset & w <= (trial.event.onset+RDK.event.duration*scr.refrate)) & ... % time of event?
                r == trial.event.RDK                    % current RDK is event RDK
            % weigh all colors by weights
            colmat_r(:,:,w) = ...
                repmat((weight(w)*RDK.RDK(r).col(1,:) + (1-weight(w))*RDK.RDK(r).col(2,:))',[1,RDK.RDK(r).num]); %Formula to compute weighted colour for every frame
            % define only for a defined proportion of dots
            colmat_r(:,t.cohidx2,w) = ...
                repmat((weight(w)*trial.event.color(1,:) + (1-weight(w))*trial.event.color(2,:))',[1,sum(t.cohidx2)]); %Formula to compute weighted colour for every frame
        else
            colmat_r(:,:,w) = ...
                repmat((weight(w)*RDK.RDK(r).col(1,:) + (1-weight(w))*RDK.RDK(r).col(2,:))',[1,RDK.RDK(r).num]); %Formula to compute weighted colour for every frame
        end
    end
    % do not interpolate alpha!
    % current implementation creates a flicker for RDKs really alternating between the two specified colors
    % and not between 'not displayed' and 'displayed'
    alphamat = colmat_r(4,:,:); % extract alpha values
    alphamat(alphamat>0)= 1; % look for non-zero alpha values and set them to 1
    %     alphamat(alphamat>=0)= 1; % set all transparency values to 1;
    colmat_r(4,:,:) = alphamat; % replace alpha values with adjusted alpha values

    % append everything; dots of different RDKs are merged into one dimension
    colmat = cat(2,colmat,colmat_r);
    lummat = cat(1,lummat, repmat(weight,RDK.RDK(r).num,1));
end

%% for movement information [across groups shown overlappingly] [they need to be comparable in parameters of size and movement etc.]
%%%% troubleshooting
num_iter = nan(r,frames.pertrial);

for r = 1:size(RDK.RDK,2)

    % write to new RDKt structure for treating them overlappingly
    RDKt = RDK.RDK(1,r); RDKt.num = sum([RDK.RDK(:,r).num]); % what number of stimuli
    RDKt.numRDKid = cell2mat(arrayfun(@(x,y) repmat(y,1,x), [RDK.RDK(:,r).num], [RDK.RDK(:,r).id], 'UniformOutput', false)); % which RDKid?


    % movdot specifies when to do dot position updates
    if RDKt.mov_freq == 0
        movdot = double([abs(diff([2 weight]))==~0]); % Writes a "1" for every first on-frame and every first off-frame
        % figure; plot(colidx(1,:)); hold on; plot([fr_cue fr_cue], [-0.1 2.1]); plot(movdot)
    else
        movdot_update = scr.refrate/RDKt.mov_freq;
        movdot = zeros(1,size(weight,2));
        movdot(fr_cue:movdot_update:end)=1; % again synced at cue onset
        movdot(fr_cue:-movdot_update:1)=1;
        % figure; plot(colidx(1,:)); hold on; plot([fr_cue fr_cue], [-0.1 2.1]); plot(movdot)
    end

    % initalize RDK positions
    if numel(RDKt.dot_size) == 1 % square dot moving?
        RDKt.rad_x = RDKt.size(1)/2-RDKt.dot_size/2; % Define the maximal x and y positions for a dot
        RDKt.rad_y = RDKt.size(2)/2-RDKt.dot_size/2;
    elseif numel(RDK.RDK(r).dot_size) == 2 % rectangle dot moving?
        RDKt.rad_x = RDKt.size(1)/2-ceil(sqrt(sum((RDKt.dot_size).^2))/2); % Define the maximal x and y positions for a rectangle (based on pythagoral distance)
        RDKt.rad_y = RDKt.size(2)/2-ceil(sqrt(sum((RDKt.dot_size).^2))/2);
    end

    if numel(RDKt.dot_size) == 1 % square dot moving?
        cross.rad = (cross.size/2+RDKt.dot_size)*cross.cutout; % Define the minimal x and y positions for a dot (cutout around fixation cross) %ToDo: is it ok to calculate with inf? cause cross.rad may be 0 and is a divisor later
    elseif numel(RDKt.dot_size) == 2 % rectangle dot moving?
        cross.rad = (cross.size/2+ceil(sqrt(sum((RDKt.dot_size).^2))/2))*cross.cutout; % Define the minimal x and y positions for a dot (cutout around fixation cross) %ToDo: is it ok to calculate with inf? cause cross.rad may be 0 and is a divisor later
    end
    dot_pos = [randi([-RDKt.rad_x RDKt.rad_x],RDKt.num,1) randi([-RDKt.rad_y RDKt.rad_y],RDKt.num,1)]'; % Generate initial dot position in a square space
    % check whether dots are within RDK area
    [dot_outside] = check_RDKarea(dot_pos,RDKt);
    [dot_inside] = check_crosscutout(dot_pos,RDKt,cross);
    [dot_overlap] = check_dotoverlap(dot_pos,RDKt);
    
    % include some fail switch
    redo_flag = 1;
    maxiter = 200;
    while redo_flag==1

        dot_reinitialize = any([dot_outside.any; dot_inside; dot_overlap],1);
        while any(dot_reinitialize) % generates new positions for only the invalid dots
            dot_pos(:,dot_reinitialize) = [randi([- RDKt.rad_x  RDKt.rad_x],sum(dot_reinitialize),1) randi([- RDKt.rad_y  RDKt.rad_y],sum(dot_reinitialize),1)]';
            [dot_outside] = check_RDKarea(dot_pos, RDKt);
            [dot_inside] = check_crosscutout(dot_pos, RDKt,cross);
            [dot_overlap] = check_dotoverlap(dot_pos,RDKt);
            dot_reinitialize = any([dot_outside.any; dot_inside; dot_overlap],1);
        end

        % setup dot position matrix
        dotmat_r = zeros(2,  RDKt.num, frames.pertrial); %setup dot matrix: xy x dot number x frames
        dot_pos_n = dot_pos;

        for f = 1:frames.pertrial %predefine dot positions for every frame

            % randomize directions first
            t.rdiridx = randi(length(RDKt.mov_dir),1,RDKt.num);
            % apply to position
            dot_pos_n = dot_pos + RDKt.mov_dir(t.rdiridx,:)'*RDKt.mov_speed*movdot(f);

            % check whether dots are within RDK area
            [dot_outside,dot_crit] = check_RDKarea(dot_pos_n,RDKt);
            % check whether dots are outside fixation cross area (if necessary)
            [dot_inside,cross_crit] = check_crosscutout(dot_pos_n,RDKt,cross);
            % check whether dots overlap
            [dot_overlap] = check_dotoverlap(dot_pos_n,RDKt);

            dot_reinitialize = any([dot_outside.any; dot_inside; dot_overlap],1);
            t_num = 0;
            while any(dot_reinitialize) & t_num <= maxiter
                t.rdiridx_n = t.rdiridx;
                for i_dot = find(dot_reinitialize)
                    t.rdiridx_n(i_dot) = randsample(setdiff(1:length(RDKt.mov_dir),t.rdiridx(i_dot)),1);
                end

                % apply to position
                dot_pos_n = dot_pos + RDKt.mov_dir(t.rdiridx_n,:)'*RDKt.mov_speed*movdot(f);
                [dot_outside,dot_crit] = check_RDKarea(dot_pos_n,RDKt);
                [dot_inside,cross_crit] = check_crosscutout(dot_pos_n,RDKt,cross);
                [dot_overlap] = check_dotoverlap(dot_pos_n,RDKt);
                dot_reinitialize = any([dot_outside.any; dot_inside; dot_overlap],1);
                %%%% troubleshooting
                t_num = t_num+1; % count the number of iterations
            end
            if t_num > maxiter
                break % redo the randomization
            end

            %%%% troubleshooting
            num_iter(r,f) = t_num;

            dot_pos = dot_pos_n;
            %write into dot matrix
            dotmat_r(:,:,f) = dot_pos;
        end

        % check whether randomization needs to be redone
        if t_num > maxiter
            redo_flag = 1;
        else
            redo_flag = 0;
        end

    end

    % add QUAD4X shift to positions
    if Propixx == 4
        for p=1:Propixx
            shift_idx = p:Propixx:frames.pertrial;
            dotmat_r(:,:,shift_idx) = [dotmat_r(1,:,shift_idx)+scr.shift(p,1); dotmat_r(2,:,shift_idx)+scr.shift(p,2)];
        end
    end
    % add RDK centershift to positions
    dotmat_r = [dotmat_r(1,:,:)+RDKt.centershift(1); dotmat_r(2,:,:)+RDKt.centershift(2)];
    % figure; scatter(dotmat_r(1,:,1),dotmat_r(2,:,1))

    % setup array to mark to wich RDK a dot belongs
    rdkidx_r = repmat(RDKt.numRDKid',1,size(dotmat_r,3));

    % define size for each dot
    dotsize_r = repmat(ceil(sqrt(sum((RDKt.dot_size).^2))/2),size(rdkidx_r));

    % append everything; dots of different RDKs are merged into one dimension
    dotmat = cat(2,dotmat,dotmat_r);
    rdkidx = cat(1,rdkidx,rdkidx_r);
    dotsize = cat(1,dotsize,dotsize_r);


    % troublehsooting:
    % plotting
    %     figure; plot(squeeze(colmat(1,:,:))')


end

%% redistribute data
colmat = reshape(colmat,[size(colmat,1), size(rdkidx,1)*Propixx, frames.flips]);
dotmat = reshape(dotmat,[2, size(rdkidx,1)*Propixx, frames.flips]);
rdkidx = reshape(rdkidx,[size(rdkidx,1)*Propixx, frames.flips]);
dotsize = reshape(dotsize,[size(dotsize,1)*Propixx, frames.flips]);
lummat = reshape(lummat,[size(lummat,1)*Propixx, frames.flips]);


varargout{1} = lummat;

%%%troubleshooting
% figure; plot(num_iter')
end

%% SUBFUNCTIONS

function [dot_outside,dot_crit] = check_RDKarea(dot_pos,RDK)
dot_crit.x = dot_pos(1,:).^2/RDK.rad_x^2;
dot_crit.y = dot_pos(2,:).^2/RDK.rad_y^2;
dot_crit.xy = dot_crit.x + dot_crit.y;

dot_outside.ellipse = dot_crit.xy > 1;
dot_outside.x = dot_crit.x > 1;
dot_outside.y = dot_crit.y > 1;
dot_outside.square = dot_outside.x | dot_outside.y;

switch RDK.shape
    case 1
        dot_outside.any = dot_outside.square;
    case 0
        dot_outside.any = dot_outside.ellipse;
end

end

function [dot_inside,cross_crit] = check_crosscutout(dot_pos,RDK,cross)
cross_crit = (dot_pos(1,:)+RDK.centershift(1)).^2/cross.rad^2 ...
    + (dot_pos(2,:)+RDK.centershift(2)).^2/cross.rad^2; % dots inside the fixation cross cutout (if cross.cutout nonzero); Note: these calculations still assume that RDK is central, centershift is added later; that's why we need to add centershift here for checking;
dot_inside = cross_crit <= 1;
end

function [dot_overlap] = check_dotoverlap(dot_pos,RDK)
% critical distance between points
if numel(RDK.dot_size) == 1 % square dot moving?
    critdist = RDK.dot_size/2;
elseif numel(RDK.dot_size) == 2 % square dot moving?
    critdist = ceil(sqrt(sum((RDK.dot_size).^2))/2);
end
dot_crit.x = abs(dot_pos(1,:)'-dot_pos(1,:)).^2/critdist^2;
dot_crit.y = abs(dot_pos(2,:)'-dot_pos(2,:)).^2/critdist^2;

% find dots too close to each other
critmat = (dot_crit.x <=1 | dot_crit.y <= 1) & tril(~eye(size(dot_crit.x,1))); % looks for any differences in matrix in lower triangular

% index the dots
[critdot1, critdot2] = ind2sub(size(critmat),find(critmat));

% find dots too close to each other
dot_overlap = false(1,size(dot_pos,2));
dot_overlap(unique([critdot1 critdot2])) = true;
end
