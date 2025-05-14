function [] = pres_feedback(responses_in,p,ps, key,RDK)
%PRES_FEEDBACK calculate and display feedback for previous block
%   Returns Percentage hit, false alarms, and reaction time

WaitSecs(0.5);
%% calculate feedback [for the all blocks]
summcon = [];
% loop across blocks
for i_block = 1:numel(responses_in)
    responses = responses_in{i_block};
    % get number of all events
    t.num_presses=sum(cellfun(@(x) sum(sum(~isnan(x))),{responses.button_presses_t}));    % number of total button presse

    % target number
    summ.targnum(i_block) = sum(cellfun(@(x) sum(x==1),{responses.eventtype}));
    summ.distrnum(i_block) = sum(cellfun(@(x) sum(x==2),{responses.eventtype}));

    summ.hits(i_block) = sum(cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.event_response_type},'UniformOutput',false)));
    summ.misses(i_block) = sum(cell2mat(cellfun(@(x) strcmpi(x,'miss'),{responses.event_response_type},'UniformOutput',false)));
    summ.CR(i_block) = sum(cell2mat(cellfun(@(x) strcmpi(x,'CR'),{responses.event_response_type},'UniformOutput',false)));
    summ.FA_proper(i_block) = sum(cell2mat(cellfun(@(x) strcmpi(x,'FA_proper'),{responses.event_response_type},'UniformOutput',false)));
    t.FA = [];
    for i_tr = 1:numel(responses)
        t.FA(i_tr) = sum(strcmp(responses(i_tr).button_presses_type,'FA'));
    end
    summ.FA(i_block) = sum(t.FA);
    summ.RT_mean(i_block) = nanmean(cell2mat(arrayfun(@(x,y) x(y), cell2mat({responses.event_response_RT}),...
        cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.event_response_type},'UniformOutput',false)),'UniformOutput',false)));
    summ.RT_std(i_block) = nanstd(cell2mat(arrayfun(@(x,y) x(y), cell2mat({responses.event_response_RT}),...
        cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.event_response_type},'UniformOutput',false)),'UniformOutput',false)));

    summ.precue_eventnum(i_block) = sum([responses.precue_eventnum]);
    summ.precue_targnum(i_block) = sum([responses.precue_eventtype]==1);
    summ.precue_distrnum(i_block) = sum([responses.precue_eventtype]==2);
    summ.precue_hits(i_block) = sum(cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.precue_event_response_type},'UniformOutput',false)));
    summ.precue_misses(i_block) = sum(cell2mat(cellfun(@(x) strcmpi(x,'miss'),{responses.precue_event_response_type},'UniformOutput',false)));
    summ.precue_error(i_block) = sum(cell2mat(cellfun(@(x) strcmpi(x,'error'),{responses.precue_event_response_type},'UniformOutput',false)));
    summ.precue_RT_mean(i_block) = nanmean(cell2mat(arrayfun(@(x,y) x(y), cell2mat({responses.precue_event_response_RT}),...
        cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.precue_event_response_type},'UniformOutput',false)),'UniformOutput',false)));
    summ.precue_RT_std(i_block) = nanstd(cell2mat(arrayfun(@(x,y) x(y), cell2mat({responses.precue_event_response_RT}),...
        cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.precue_event_response_type},'UniformOutput',false)),'UniformOutput',false)));


    
    % behavioral effects separately for experimental conditions
    for i_con = 1:numel(p.stim.condition)
        % extract indices
        t.idx_con = repmat([responses.cue]==p.stim.condition(i_con),size([responses.eventtype2],1),1);
        % extract responses for cue
        t.responses = cellfun(@(x) x',{responses.event_response_type},'UniformOutput',false);
        t.responses = [t.responses{:}];

        % extract response times
        t.response_RTs = reshape([responses.event_response_RT],size(t.responses,1),[]);
        for i_evtype = 1:2 % primed or non primed event
            % how many (different) targets?
            summcon.targnum(i_con,i_evtype,i_block) = ...
                sum([responses.eventtype2]==i_evtype & t.idx_con,"all");
            % index primed or non-primed event
            t.idx_evtype_t = [responses.eventtype2]==i_evtype;
            summcon.hits(i_con,i_evtype,i_block) = sum(strcmp(t.responses,'hit')&t.idx_con&t.idx_evtype_t,'all');
            summcon.misses(i_con,i_evtype,i_block) = sum(strcmp(t.responses,'miss')&t.idx_con&t.idx_evtype_t,'all');

            summcon.RT{i_con,i_evtype,i_block} = t.response_RTs(strcmp(t.responses,'hit')&t.idx_con&t.idx_evtype_t)';
        end
        % how many distractors?
        t.idx_evtype_d = [responses.eventtype2]==3;
        summcon.distrnum(i_con,i_block) = sum(t.idx_evtype_d & t.idx_con,"all");
        % extract FA_proper and CR
        summcon.CR(i_con,i_block) = sum(strcmp(t.responses(t.idx_evtype_d & t.idx_con),'CR'),"all");
        summcon.FA_proper(i_con,i_block) = sum(strcmp(t.responses(t.idx_evtype_d & t.idx_con),'FA_proper'),"all");

        % extract FAs
        t.FA = [];
        t.idxFA = find([responses.cue]==p.stim.condition(i_con));
        for i_tr = 1:numel(t.idxFA)
            t.FA(i_tr) = sum(strcmp(responses(t.idxFA(i_tr)).button_presses_type,'FA'));
        end
        summcon.FA(i_con,i_block) = sum(t.FA);
    end
end

%% data needs to be aggregated
% show last block + all so far run blocks
% show separate for the different target types (primed + non-primed) for participants
% show separately for conditions outside as well
try
dat2disp.all.Hitrate = [summ.hits(end)/summ.targnum(end) sum(summ.hits)/sum(summ.targnum)]*100; %[last block, overall]
dat2disp.all.FArate = [summ.FA_proper(end)/summ.distrnum(end) sum(summ.FA_proper)/sum(summ.distrnum)]*100; %[last block, overall]
dat2disp.all.targets = [summ.targnum(end) sum(summ.targnum)]; %[last block, overall]
dat2disp.all.distractors = [summ.distrnum(end) sum(summ.distrnum)]; %[last block, overall]
dat2disp.all.RT_mean = [mean([summcon.RT{:,:,end}],"all") mean([summcon.RT{:}],"all")];
dat2disp.all.RT_std = [std([summcon.RT{:,:,end}],0,"all") std([summcon.RT{:}],0,"all")];

dat2disp.precue.Hitrate = [summ.precue_hits(end)/summ.precue_targnum(end) sum(summ.precue_hits)/sum(summ.precue_targnum)]*100; %[last block, overall]
dat2disp.precue.FArate = [summ.precue_error(end)/summ.precue_distrnum(end) sum(summ.precue_error)/sum(summ.precue_distrnum)]*100; %[last block, overall]
dat2disp.precue.targets = [summ.precue_targnum(end) sum(summ.precue_targnum)]; %[last block, overall]
dat2disp.precue.distractors = [summ.precue_distrnum(end) sum(summ.precue_distrnum)]; %[last block, overall]
dat2disp.precue.RT_mean = [summ.precue_RT_mean(end) mean(summ.precue_RT_mean)];
dat2disp.precue.RT_std = [summ.precue_RT_std(end) mean(summ.precue_RT_std)];

dat2disp.byprime.Hitrate = [(sum(summcon.hits(:,:,end))./sum(summcon.targnum(:,:,end)))' ... % sum of hits for [primed non-primed] event
    (sum(summcon.hits,[1,3])./sum(summcon.targnum,[1,3]))']*100; %[last block primed, overall primed; last block non-primed, overall non-primed]
dat2disp.byprime.targets = [sum(summcon.targnum(:,:,end))' sum(summcon.targnum,[1,3])']; %[last block, overall]
dat2disp.byprime.RT_mean = [mean([summcon.RT{:,1,end}],"all") mean([summcon.RT{:,1,:}],"all"); ...
    mean([summcon.RT{:,2,end}],"all") mean([summcon.RT{:,2,:}],"all")];
dat2disp.byprime.RT_std = [std([summcon.RT{:,1,end}],1,"all") std([summcon.RT{:,1,:}],1,"all"); ...
    std([summcon.RT{:,2,end}],1,"all") std([summcon.RT{:,2,:}],1,"all")];

dat2disp.bycon.last.Hitrate = (summcon.hits(:,:,end)./summcon.targnum(:,:,end))*100; % separate by prime and non-prime
dat2disp.bycon.last.FArate = (summcon.FA_proper(:,end)./summcon.distrnum(:,end))*100;
dat2disp.bycon.last.targets = summcon.targnum(:,:,end); % separate by prime and non-prime
dat2disp.bycon.last.distractors = summcon.distrnum(:,end); % separate by prime and non-primed
[dat2disp.bycon.last.RT_mean, dat2disp.bycon.last.RT_std] = deal([]);
for i_con = 1:numel(p.stim.condition)
    for i_prime = 1:2
        dat2disp.bycon.last.RT_mean(i_con,i_prime) = mean([summcon.RT{i_con,i_prime,end}],"all");
        dat2disp.bycon.last.RT_std(i_con,i_prime) = std([summcon.RT{i_con,i_prime,end}],1,"all");
    end
end

dat2disp.bycon.all.Hitrate = (sum(summcon.hits,3)./sum(summcon.targnum,3))*100; % separate by prime and non-prime
dat2disp.bycon.all.FArate = (sum(summcon.FA_proper,2)./sum(summcon.distrnum,2))*100;
dat2disp.bycon.all.targets = sum(summcon.targnum,3); % separate by prime and non-prime
dat2disp.bycon.all.distractors = sum(summcon.distrnum,2); % separate by prime and non-primed
[dat2disp.bycon.all.RT_mean, dat2disp.bycon.all.RT_std] = deal([]);
for i_con = 1:numel(p.stim.condition)
    for i_prime = 1:2
        dat2disp.bycon.all.RT_mean(i_con,i_prime) = mean([summcon.RT{i_con,i_prime,:}],"all");
        dat2disp.bycon.all.RT_std(i_con,i_prime) = std([summcon.RT{i_con,i_prime,:}],0,"all");
    end
end
catch
    fprintf('error')
end



%% presentation of results
% KbQueueCreate(ps.RespDev, keysOfInterest) %
% KbQueueStart(ps.RespDev);

% pixels for shift into 4 quadrants
quadshift = [p.scr_res(1)*(1/4) p.scr_res(2)*(1/4); p.scr_res(1)*(3/4) p.scr_res(2)*(1/4); ...
    p.scr_res(1)*(1/4) p.scr_res(2)*(3/4); p.scr_res(1)*(3/4) p.scr_res(2)*(3/4)];

% output to screen
fprintf('\n###block %02d###',numel(responses_in))
% fprintf(['\nPre-Cue:\n Hitrate: %06.2f(%06.2f)%%; targs: %03d(%03d); FA_proper: %06.2f(%06.2f)%%; distr: %03d(%03d); RT_M: %3.0f(%3.0f), RT_Std: %3.0f(%3.0f) ms'], ...
%     dat2disp.precue.Hitrate, dat2disp.precue.targets, dat2disp.precue.FArate, dat2disp.precue.distractors,dat2disp.precue.RT_mean, dat2disp.precue.RT_std)
fprintf(['\nMAIN | all collapsed:\n Hitrate: %06.2f(%06.2f)%%; targs: %03d(%03d); FA_proper: %06.2f(%06.2f)%%; distr: %03d(%03d); RT_M: %3.0f(%3.0f), RT_Std: %3.0f(%3.0f) ms'], ...
    dat2disp.all.Hitrate, dat2disp.all.targets, dat2disp.all.FArate, dat2disp.all.distractors,dat2disp.all.RT_mean, dat2disp.all.RT_std)
fprintf('\nMAIN by prime | color collapsed:')
fprintf(['\n primed:    Hitrate: %06.2f(%06.2f)%%; targs: %03d(%03d); RT_M: %3.0f(%3.0f), RT_Std: %3.0f(%3.0f) ms'], ...
    dat2disp.byprime.Hitrate(1,:), dat2disp.byprime.targets(1,:), dat2disp.byprime.RT_mean(1,:), dat2disp.byprime.RT_std(1,:))
fprintf(['\n nonprimed: Hitrate: %06.2f(%06.2f)%%; targs: %03d(%03d); RT_M: %3.0f(%3.0f), RT_Std: %3.0f(%3.0f) ms'], ...
    dat2disp.byprime.Hitrate(2,:), dat2disp.byprime.targets(2,:), dat2disp.byprime.RT_mean(2,:), dat2disp.byprime.RT_std(2,:))
for i_con = 1:numel(p.stim.condition)
    fprintf('\nMAIN by prime | attend [%s %s]:', RDK.RDK(p.stim.RDK2attend(i_con,:)).col_label)
    fprintf(['\n primed:    Hitrate: %06.2f(%06.2f)%%; targs: %03d(%03d); RT_M: %3.0f(%3.0f), RT_Std: %3.0f(%3.0f) ms'], ...
        dat2disp.bycon.last.Hitrate(i_con,1), dat2disp.bycon.all.Hitrate(i_con,1), ...
        dat2disp.bycon.last.targets(i_con,1), dat2disp.bycon.all.targets(i_con,1), ...
        dat2disp.bycon.last.RT_mean(i_con,1), dat2disp.bycon.all.RT_mean(i_con,1), ...
        dat2disp.bycon.last.RT_std(i_con,1), dat2disp.bycon.all.RT_std(i_con,1))
    fprintf(['\n nonprimed: Hitrate: %06.2f(%06.2f)%%; targs: %03d(%03d); RT_M: %3.0f(%3.0f), RT_Std: %3.0f(%3.0f) ms'], ...
        dat2disp.bycon.last.Hitrate(i_con,2), dat2disp.bycon.all.Hitrate(i_con,2), ...
        dat2disp.bycon.last.targets(i_con,2), dat2disp.bycon.all.targets(i_con,2), ...
        dat2disp.bycon.last.RT_mean(i_con,2), dat2disp.bycon.all.RT_mean(i_con,2), ...
        dat2disp.bycon.last.RT_std(i_con,2), dat2disp.bycon.all.RT_std(i_con,2))
end
fprintf('\n###\nWeiter mit q...\n')




% draw text and stimuli (before shifting to the quadrants)
% offscreen window
[ps.offwin,ps.offrect]=Screen('OpenOffscreenWindow',p.scr_num, [0 0 0 0], [0 0 p.scr_res(1)/2 p.scr_res(2)/2], [], [], []);
% get center of offscreen window
[ps.xCenter_off, ps.yCenter_off] = RectCenter(ps.offrect);

% text2present=                   [...                % text for feedback
%     'P A U S E'...
%     sprintf('\n\nmain task | letzter Block (alle)')...
%     sprintf('\nHitrate   = %06.2f(%06.2f)%%,       Zielreize = %03d(%03d)',dat2disp.all.Hitrate, dat2disp.all.targets)...
%     sprintf('\nFA_rate = %06.2f(%06.2f)%%, Distraktoren = %03d(%03d)',dat2disp.all.FArate, dat2disp.all.distractors)...
%     sprintf('\nReaktionszeit: M = %1.0f(%1.0f)ms, Std = %1.0f(%1.0f)ms',dat2disp.all.RT_mean, dat2disp.all.RT_std)...
%     sprintf('\n\npre-cue task | letzter Block (alle)')...
%     sprintf('\nHitrate   = %06.2f(%06.2f)%%,       Zielreize = %03d(%03d)', dat2disp.precue.Hitrate, dat2disp.precue.targets)...
%     sprintf('\nFA_rate = %06.2f(%06.2f)%%, Distraktoren = %03d(%03d)',dat2disp.precue.FArate, dat2disp.precue.distractors)...
%     sprintf('\nReaktionszeit: M = %1.0f(%1.0f)ms, Std = %1.0f(%1.0f)ms',dat2disp.precue.RT_mean, dat2disp.precue.RT_std)...
%     ];

% text2present=                   [...                % text for feedback
%     'P A U S E'...
%     sprintf('\n\nmain task | letzter Block (alle)')...
%     sprintf('\nHitrate   = %06.2f(%06.2f)%%',dat2disp.all.Hitrate)...
%     sprintf('\nFA_rate = %06.2f(%06.2f)%%',dat2disp.all.FArate)...
%     sprintf('\n\n\npre-cue task | letzter Block (alle)')...
%     sprintf('\nHitrate   = %06.2f(%06.2f)%%', dat2disp.precue.Hitrate)...
%     sprintf('\nFA_rate = %06.2f(%06.2f)%%',dat2disp.precue.FArate)...
%     ];

% no pre-cue task and only reaction times?
text2present=                   [...                % text for feedback
    'P A U S E'...
    sprintf('\nReaktionszeit: M = %1.0f(%1.0f)ms, Std = %1.0f(%1.0f)ms',dat2disp.all.RT_mean, dat2disp.all.RT_std)...
    sprintf('\nauf %1.0f(%1.0f) Zielreize reagiert von %1.0f(%1.0f) möglichen', ...
    dat2disp.all.Hitrate./100.*dat2disp.all.targets, dat2disp.all.targets)
    ];


% draw text
Screen('TextSize', ps.offwin, 18);
% DrawFormattedText(tx.instruct, INST.text{1}, tx.xCenter_off, p.scr_res(1)/2 * 0.1, p.stim_color);
DrawFormattedText(ps.offwin, text2present, 'center', p.scr_res(1)/2 * 0.2, p.crs.color);

for i_quad = 1:4 % shifst to quadrants
    newpos_stim(:,i_quad) = ...
        CenterRectOnPointd(ps.offrect,quadshift(i_quad,1),quadshift(i_quad,2))';
end
 
% [key.pressed, key.firstPress]=KbQueueCheck;
key.rkey=key.SECRET;
[key.keyisdown,key.secs,key.keycode] = KbCheck; 
while ~(key.keycode(key.rkey)==1)                       % continuously present feedback (wait for q)
    [key.keyisdown,key.secs,key.keycode] = KbCheck;
    Screen('DrawTextures', ps.window, repmat(ps.offwin,1,4),[], newpos_stim, [], [], [], []);
    Screen('Flip', ps.window, 0);                   % flip screen
end
try Screen('Close', ps.offwin); end


end

