function [] = pres_feedback(responses,p,ps, key,RDK)
%PRES_FEEDBACK calculate and display feedback for previous block
%   Returns Percentage hit, false alarms, and reaction time

WaitSecs(0.5);
%% calculate feedback
% get number of all events
t.num_presses=sum(cellfun(@(x) sum(sum(~isnan(x))),{responses.button_presses_t}));    % number of total button presse

summ.hits = sum(cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.event_response_type},'UniformOutput',false)));
summ.misses = sum(cell2mat(cellfun(@(x) strcmpi(x,'miss'),{responses.event_response_type},'UniformOutput',false)));
summ.error = sum(cell2mat(cellfun(@(x) strcmpi(x,'error'),{responses.event_response_type},'UniformOutput',false)));
summ.RT_mean = nanmean(cell2mat(arrayfun(@(x,y) x(y), cell2mat({responses.event_response_RT}),...
    cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.event_response_type},'UniformOutput',false)),'UniformOutput',false)));
summ.RT_std = nanstd(cell2mat(arrayfun(@(x,y) x(y), cell2mat({responses.event_response_RT}),...
    cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.event_response_type},'UniformOutput',false)),'UniformOutput',false)));
summ.totalevents = numel({responses.event_response_type});

summ.precue_eventnum = sum([responses.precue_eventnum]);
summ.precue_hits = sum(cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.precue_event_response_type},'UniformOutput',false)));
summ.precue_misses = sum(cell2mat(cellfun(@(x) strcmpi(x,'miss'),{responses.precue_event_response_type},'UniformOutput',false)));
summ.precue_error = sum(cell2mat(cellfun(@(x) strcmpi(x,'error'),{responses.precue_event_response_type},'UniformOutput',false)));
summ.precue_RT_mean = nanmean(cell2mat(arrayfun(@(x,y) x(y), cell2mat({responses.precue_event_response_RT}),...
    cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.precue_event_response_type},'UniformOutput',false)),'UniformOutput',false)));
summ.precue_RT_std = nanstd(cell2mat(arrayfun(@(x,y) x(y), cell2mat({responses.precue_event_response_RT}),...
    cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.precue_event_response_type},'UniformOutput',false)),'UniformOutput',false)));

% behavioral effects separately for experimental conditions
for i_con = 1:6
    summ.hits(1+i_con) = sum(strcmp({responses.event_response_type},'hit') & [responses.condition]==i_con);
    summ.misses(1+i_con) = sum(strcmp({responses.event_response_type},'miss') & [responses.condition]==i_con);
    summ.error(1+i_con) = sum(strcmp({responses.event_response_type},'error') & [responses.condition]==i_con);

    summ.RT_mean(1+i_con) = mean([responses( ...
        strcmp({responses.event_response_type},'hit') & [responses.condition]==i_con) ...
        .event_response_RT],"omitnan");
    
    summ.RT_std(1+i_con) = std([responses( ...
        strcmp({responses.event_response_type},'hit') & [responses.condition]==i_con) ...
        .event_response_RT],"omitnan");

    summ.totalevents(1+i_con) = sum([responses.condition]==i_con);
end

% pixels for shift into 4 quadrants
quadshift = [p.scr_res(1)*(1/4) p.scr_res(2)*(1/4); p.scr_res(1)*(3/4) p.scr_res(2)*(1/4); ...
    p.scr_res(1)*(1/4) p.scr_res(2)*(3/4); p.scr_res(1)*(3/4) p.scr_res(2)*(3/4)];

%% set up display of results
% KbQueueCreate(ps.RespDev, keysOfInterest) %
% KbQueueStart(ps.RespDev);

% output to screen
t.textin = {...
    'All:           ';...
    'cue_left_3t3d  ';
    'cue_right_3t3d ';
    'cue_left_6t    ';
    'cue_right_6t   ';
    'cue_left_3t    ';
    'cue_right_3t   ';};
for i_con = 1:numel(summ.hits)
    fprintf(1,...
        ['\n%s Hitrate: %06.2f; Hits: %02.0f; Misses: %02.0f; Error: %02.0f; RT: M: %3.0f, Std: %3.0f ms'],...
        t.textin{i_con}, summ.hits(i_con)/summ.totalevents(i_con)*100, summ.hits(i_con), summ.misses(i_con), ...
        summ.error(i_con), summ.RT_mean(i_con), summ.RT_std(i_con))
end
% pre-cue events
fprintf('\npre-cue events: Hits: %02.0f; Misses: %02.0f; Error: %02.0f; RT: M: %3.0f, Std: %3.0f ms', ...
    summ.precue_hits, summ.precue_misses, summ.precue_error, summ.precue_RT_mean, summ.precue_RT_std)
fprintf('\nExperiment freigeben mit "q".')


% display 1
% draw text and stimuli (before shifting to the quadrants) for feedback
% offscreen window
[ps.offwin,ps.offrect]=Screen('OpenOffscreenWindow',p.scr_num, [0 0 0 0], [0 0 p.scr_res(1)/2 p.scr_res(2)/2], [], [], []);
% get center of offscreen window
[ps.xCenter_off, ps.yCenter_off] = RectCenter(ps.offrect);

text2present=                   [...                % text for feedback
    'P A U S E'...
    sprintf('\n\nHits = %1.0f; Hitrate = %1.2f%%',summ.hits(1), summ.hits(1)/summ.totalevents(1)*100)...
    sprintf('\n\nMisses = %1.0f; Error = %1.0f',summ.misses(1), summ.error(1))...
    sprintf('\n\nReaktionszeit: M = %1.0fms, Std = %1.0fms',summ.RT_mean(1), summ.RT_std(1))...
    sprintf('\n\n\npre-cue: Hits = %1.0f; Misses = %1.0f; Error = %1.0f',	summ.precue_hits, summ.precue_misses, summ.precue_error)...
    ];

% draw text
Screen('TextSize', ps.offwin, 18);
% DrawFormattedText(tx.instruct, INST.text{1}, tx.xCenter_off, p.scr_res(1)/2 * 0.1, p.stim_color);
DrawFormattedText(ps.offwin, text2present, 'center', p.scr_res(1)/2 * 0.2, p.crs.color);

for i_quad = 1:4 % shifst to quadrants
    newpos_stim(:,i_quad) = ...
        CenterRectOnPointd(ps.offrect,quadshift(i_quad,1),quadshift(i_quad,2))';
end

% display 2
% draw text and stimuli (before shifting to the quadrants) for feedback
% offscreen window
[ps.offwin2,ps.offrect2]=Screen('OpenOffscreenWindow',p.scr_num, [0 0 0 0], [0 0 p.scr_res(1)/2 p.scr_res(2)/2], [], [], []);
% get center of offscreen window
[ps.xCenter_off, ps.yCenter_off] = RectCenter(ps.offrect2);

text2present2=                   [...                % text for feedback
    'P A U S E'...
    sprintf('\n\nHits = %1.0f; Hitrate = %1.2f%%',summ.hits(1), summ.hits(1)/summ.totalevents(1)*100)...
    sprintf('\n\nMisses = %1.0f; Error = %1.0f',summ.misses(1), summ.error(1))...
    sprintf('\n\nReaktionszeit: M = %1.0fms, Std = %1.0fms',summ.RT_mean(1), summ.RT_std(1))...
    sprintf('\n\n\npre-cue: Hits = %1.0f; Misses = %1.0f; Error = %1.0f',	summ.precue_hits, summ.precue_misses, summ.precue_error)...
    '\n\nLeertaste drücken um fortzufahren'
    ];

% draw text
Screen('TextSize', ps.offwin2, 18);
% DrawFormattedText(tx.instruct, INST.text{1}, tx.xCenter_off, p.scr_res(1)/2 * 0.1, p.stim_color);
DrawFormattedText(ps.offwin2, text2present2, 'center', p.scr_res(1)/2 * 0.2, p.crs.color);

 
% [key.pressed, key.firstPress]=KbQueueCheck;
% first check for button press of Q
key.rkey=key.SECRET;
[key.keyisdown,key.secs,key.keycode] = KbCheck; 
while ~(key.keycode(key.rkey)==1)                       % continuously present feedback (wait for q)
    [key.keyisdown,key.secs,key.keycode] = KbCheck;
    Screen('DrawTextures', ps.window, repmat(ps.offwin,1,4),[], newpos_stim, [], [], [], []);
    Screen('Flip', ps.window, 0);                   % flip screen
end

% display instruction outside
fprintf('..Versuchsperson startet mit "Leertaste".\n')

% then check for button press of space
key.rkey=key.SPACE;
[key.keyisdown,key.secs,key.keycode] = KbCheck; 
while ~(key.keycode(key.rkey)==1)                       % continuously present feedback (wait for q)
    [key.keyisdown,key.secs,key.keycode] = KbCheck;
    Screen('DrawTextures', ps.window, repmat(ps.offwin2,1,4),[], newpos_stim, [], [], [], []);
    Screen('Flip', ps.window, 0);                   % flip screen
end

try Screen('Close', ps.offwin); end


end

