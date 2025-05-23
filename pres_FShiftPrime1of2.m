function [timing,key,resp] = pres_FShiftPrime1of2(p, ps, key, RDK, conmat, blocknum, flag_training)
% presents experiment SSVEP_FShiftPrime1of2
%   p               = parameters
%   ps              = screen parameters
%   RDK             = RDK parameters
%   blocknum        = number of block
%   flag_training   = flag for trainig (1) or experiment (0)

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
RDKin.crs.size = RDKin.crs.dims;

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
crttime = GetSecs;

% loop across trials
for i_tr = 1:numel(trialindex)
    fprintf('%1.0f',mod(i_tr,10))
    inittime=GetSecs;
    %% initialize trial structure, RDK, cross, logs
    RDKin.trial = struct('duration',conmat.trials(trialindex(i_tr)).post_cue_times+conmat.trials(trialindex(i_tr)).pre_cue_times,...
        'frames',conmat.trials(trialindex(i_tr)).post_cue_frames+conmat.trials(trialindex(i_tr)).pre_cue_frames,...
        'cue',conmat.trials(trialindex(i_tr)).pre_cue_frames+1);
    RDKin.trial.event = struct('onset',conmat.trials(trialindex(i_tr)).event_onset_frames,...
        'direction',conmat.trials(trialindex(i_tr)).eventdirection,'RDK',conmat.trials(trialindex(i_tr)).eventRDK);
    [colmat,dotmat,dotsize,rdkidx,frames] = RDK_init_FShiftPrime1of2(RDKin.scr,RDKin.Propixx,RDKin.RDK,RDKin.trial,RDKin.crs);
    
    % initialize fixation cross
    % which cross?
    crossmat = ones(1,size(colmat,3)); % by default show normal fixation cross
    if conmat.trials(trialindex(i_tr)).precue_eventnum == 1
        % depict fixation cross event for relevant flips (not frames)
        crossmat(conmat.trials(trialindex(i_tr)).precue_event_onset_frames/p.scr_imgmultipl-1+(1:p.scr_refrate/p.scr_imgmultipl*p.stim.precue_event.length)) = ...
            conmat.trials(trialindex(i_tr)).precue_eventid+1;
    end
    % color
    colmat_cr = repmat(p.crs.color' ,[1 1 size(colmat,3)]);
    % from cue onwards: color of to be attended RDK
    colmat_cr(:,:,(conmat.trials(trialindex(i_tr)).pre_cue_frames/p.scr_imgmultipl)+1:end)=...
        repmat(RDK.RDK(conmat.trials(trialindex(i_tr)).cue).col(1,:)',[1,1,conmat.trials(trialindex(i_tr)).post_cue_frames/p.scr_imgmultipl]);
    
    % preallocate timing
    timing(i_tr) = struct('VBLTimestamp',NaN(1,frames.flips),'StimulusOnsetTime',NaN(1,frames.flips),...
        'FlipTimestamp',NaN(1,frames.flips),'Missed',NaN(1,frames.flips));
    
    %% set up responses
    %setup key presses
    key.presses{i_tr}=nan(size(colmat,3),sum(key.keymap));
    key.presses_t{i_tr}=nan(size(colmat,3),sum(key.keymap));
    
    resp(i_tr).trialnumber              = trialindex(i_tr);
    resp(i_tr).blocknumber              = conmat.trials(trialindex(i_tr)).block;
    resp(i_tr).cue                      = conmat.trials(trialindex(i_tr)).cue; % attended RDK
    resp(i_tr).RDK_attended             = conmat.trials(trialindex(i_tr)).RDK2attend;
    resp(i_tr).RDK_unattended           = setdiff(1:numel(RDK.RDK),resp(i_tr).RDK_attended);
    resp(i_tr).col_attended             = {RDK.RDK(conmat.trials(trialindex(i_tr)).RDK2attend).col_label};
    resp(i_tr).freq_attended            = [RDK.RDK(conmat.trials(trialindex(i_tr)).RDK2attend).freq];
    resp(i_tr).cue_onset_fr             = conmat.trials(trialindex(i_tr)).pre_cue_frames + 1;
    resp(i_tr).cue_onset_t_est          = (conmat.trials(trialindex(i_tr)).pre_cue_frames + 1)/p.scr_refrate*1000;
    resp(i_tr).cue_onset_t_meas         = nan; % measured onset time for cue
    resp(i_tr).pre_cue_frames           = conmat.trials(trialindex(i_tr)).pre_cue_frames;
    resp(i_tr).pre_cue_times            = conmat.trials(trialindex(i_tr)).pre_cue_times;
    resp(i_tr).post_cue_times           = conmat.trials(trialindex(i_tr)).post_cue_frames;
    resp(i_tr).post_cue_frames          = conmat.trials(trialindex(i_tr)).post_cue_times;
    resp(i_tr).eventnum                 = conmat.trials(trialindex(i_tr)).eventnum;
    resp(i_tr).eventtype                = conmat.trials(trialindex(i_tr)).eventtype; % 1 = target; 2 = distractor
    resp(i_tr).eventtype2               = conmat.trials(trialindex(i_tr)).eventtype2; % 1 = prime target; 2= non-prime target; 3 = distractor
    resp(i_tr).eventRDK                 = conmat.trials(trialindex(i_tr)).eventRDK;
    resp(i_tr).eventdirection           = conmat.trials(trialindex(i_tr)).eventdirection;
    resp(i_tr).event_onset_frames       = conmat.trials(trialindex(i_tr)).event_onset_frames;
    resp(i_tr).event_onset_times        = conmat.trials(trialindex(i_tr)).event_onset_times;
    resp(i_tr).precue_eventnum          = conmat.trials(trialindex(i_tr)).precue_eventnum;
    resp(i_tr).precue_eventid           = conmat.trials(trialindex(i_tr)).precue_eventid;
    resp(i_tr).precue_eventtype         = conmat.trials(trialindex(i_tr)).precue_eventtype;
    resp(i_tr).precue_event_onset_fr    = conmat.trials(trialindex(i_tr)).precue_event_onset_frames;
    resp(i_tr).precue_event_onset_t_est = conmat.trials(trialindex(i_tr)).precue_event_onset_times;
    
    %% set up datapixx trigger vector
    % prepare datapixx scheduler
    % trigger signal needs to encompass first frame after stimulation has ended to send trial stop trigger
    if size(colmat_cr,3)*p.scr_imgmultipl>...
            conmat.trials(trialindex(i_tr)).pre_cue_frames+conmat.trials(trialindex(i_tr)).post_cue_frames+1
        doutWave = zeros(1,size(colmat_cr,3)*p.scr_imgmultipl); % each entry corresponds to a trigger
    else
        doutWave = zeros(1,size(colmat_cr,3)*p.scr_imgmultipl+p.scr_imgmultipl); % each entry corresponds to a trigger
    end
    
    % write trigger numbers into doutwave
    % trial start trigger
    doutWave(1) = p.trig.tr_start;
    % trial end trigger (presented at first frame after stimulation)
    doutWave(size(colmat_cr,3)*p.scr_imgmultipl+1) = p.trig.tr_stop;
    % condition trigger
    resp(i_tr).triggernum = ...
        p.trig.tr_cue_type(resp(i_tr).cue); % cue to RDK1 or RDK2?
    try resp(i_tr).triggernum = resp(i_tr).triggernum + ...
            p.trig.type(1,resp(i_tr).eventnum); % number of events?
    end
    doutWave(resp(i_tr).cue_onset_fr) = resp(i_tr).triggernum;
    
    % event trigger
    for i_ev = 1:resp(i_tr).eventnum
        t.trigger = p.trig.event_type([resp(i_tr).RDK_attended' resp(i_tr).RDK_unattended] ...
            == resp(i_tr).eventRDK(i_ev));
        doutWave(resp(i_tr).event_onset_frames(i_ev)) = t.trigger;
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
        % troubleshooting
        %fprintf('\nbefore wait = %1.3f; after wait = %1.3f\n',crttime2 - crttime, GetSecs - crttime)
    end
    
    % draw fixation cross
    Screen('DrawLines', ps.window, p.crs.lines{1}, p.crs.width, p.crs.color, ps.center, 0);
    
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
    
    %% loop across trials   
    for i_fl = 1:frames.flips
        %% Drawing
        % RDK
        Screen('DrawDots', ps.window, dotmat(:,:,i_fl), dotsize(:,i_fl), colmat(:,:,i_fl), ps.center, 0, 0);
        % fixation cross
        Screen('DrawLines', ps.window, p.crs.lines{crossmat(i_fl)}, p.crs.width, colmat_cr(:,1,i_fl)', ps.center, 0);
        
        %% start trigger schedule and start listening to response device
        if i_fl == 1 % send the trigger with the start of the 1st flip
            Datapixx('RegWrVideoSync');
        end
        
        % Flip
        [timing(i_tr).VBLTimestamp(i_fl), timing(i_tr).StimulusOnsetTime(i_fl), timing(i_tr).FlipTimestamp(i_fl), timing(i_tr).Missed(i_fl)] = Screen('Flip', ps.window, 0);
        
        % send trigger/save timing/ reset timing
        if i_fl == 1
            % start trigger
            starttime=GetSecs;
            KbEventFlush(ps.RespDev); % flush keyboard
        end
        
        %% check for button presses
        [key.pressed, key.firstPress]=KbQueueCheck(ps.RespDev);
        key.presses{i_tr}(i_fl,:)=key.firstPress(key.keymap)>1;
        key.presses_t{i_tr}(i_fl,:)=(key.firstPress(key.keymap)-starttime).*key.presses{i_tr}(i_fl,:);
        if any(key.firstPress(key.keymap)>1)
            lptwrite(1,find(key.firstPress(key.keymap),1,'first'),500);
        end
    end
    %% ITI
    % draw fixation cross again
    Screen('DrawLines', ps.window, p.crs.lines{1}, p.crs.width, p.crs.color, ps.center, 0);
    
    % flip
    Screen('Flip', ps.window, 0);
    
    % get time
    crttime = GetSecs;
    
    % add waiting period to check for late button presses
    ttt=WaitSecs(p.targ_respwin(2)/1000-p.stim.event.min_offset-p.stim.event.length);
    
    % check for button presses
    [key.pressed, key.firstPress]=KbQueueCheck(ps.RespDev);
    key.presses{i_tr}(i_fl+1,:)=key.firstPress(key.keymap)>1;
    key.presses_t{i_tr}(i_fl+1,:)=(key.firstPress(key.keymap)-starttime).*key.presses{i_tr}(i_fl+1,:);
    
    
    
    
    %%%%%%%%
    % do behavioral calculation
    key.presses{i_tr}(1,:)=[];
    key.presses_t{i_tr}(1,:)=[];
    
    % get frame and timing of button press onsets
    resp(i_tr).button_presses_fr=nan(max(sum(key.presses{i_tr})),size(key.presses{i_tr},2));
    resp(i_tr).button_presses_t=resp(i_tr).button_presses_fr;
    resp(i_tr).button_presses_t_rel2cue=resp(i_tr).button_presses_fr;
    for i_bt = 1:size(key.presses{i_tr},2)
        try
            resp(i_tr).button_presses_fr(1:sum(key.presses{i_tr}(:,i_bt)),i_bt)=...
                find(key.presses{i_tr}(:,i_bt));
            resp(i_tr).button_presses_t(1:sum(key.presses{i_tr}(:,i_bt)),i_bt)=...
                key.presses_t{i_tr}(find(key.presses{i_tr}(:,i_bt)),i_bt)*1000; % in ms
            % relative to cue
            resp(i_tr).button_presses_t_rel2cue(1:sum(key.presses{i_tr}(:,i_bt)),i_bt)=...
                key.presses_t{i_tr}(find(key.presses{i_tr}(:,i_bt)),i_bt)*1000 -...
                resp(i_tr).pre_cue_times*1000; % in ms
        catch
            resp(i_tr).button_presses_fr(:,i_bt)=nan;
        end
    end

    % all relevant presses 
    t.presses = resp(i_tr).button_presses_t(:,key.keymap_ind==key.SPACE);

    % check for pre-cue events
    % check for hits {'hit','miss','error'}
    resp(i_tr).precue_event_response_type = {}; %{'hit','miss','error','CR'}
    resp(i_tr).precue_event_response_RT = []; % reaction time in ms (after event or closest to other event)

    % first define response windows
    if any(~isnan(resp(i_tr).precue_eventtype))
        t.respwin = (resp(i_tr).precue_event_onset_t_est+(p.targ_respwin/1000))*1000;
        % check for button presses in time window
        t.idx_pre = t.presses>=t.respwin(:,1) & t.presses<=t.respwin(:,2);
    else
        t.idx_pre = false(size(t.presses));
    end

       
    % not in any response window? --> miss; no reaction time
    % % troublshooting
    % if resp(i_tr).precue_eventtype == 2
    %     t.t = 1;
    % end

    if resp(i_tr).precue_eventnum < 1
        resp(i_tr).precue_event_response_type = 'noevent';
        % no reaction time
        resp(i_tr).precue_event_response_RT = nan;
    elseif ~any(t.idx_pre) & resp(i_tr).precue_eventnum == 1 & resp(i_tr).precue_eventtype == 1
        resp(i_tr).precue_event_response_type = 'miss';
        % no reaction time
        resp(i_tr).precue_event_response_RT = nan;
    elseif ~any(t.idx_pre) & resp(i_tr).precue_eventnum == 1 & resp(i_tr).precue_eventtype == 2
        resp(i_tr).precue_event_response_type = 'CR';
        % no reaction time
        resp(i_tr).precue_event_response_RT = nan;
    else
        % check for first respons (in case of two presses: should be rare but might happen)
        [~, t.idx2] = min(t.presses(t.idx_pre));
        t.idx3 = false(size(t.presses)); t.idx3(t.idx_pre(t.idx2))=true;

        if find(any(t.idx3,1)) == resp(i_tr).precue_eventtype % check for button presses of same for same trials and diff for diff trials
            resp(i_tr).precue_event_response_type = 'hit';
        else
            resp(i_tr).precue_event_response_type = 'error';
        end

        % reaction time?
        resp(i_tr).precue_event_response_RT = t.presses(t.idx3)- ...
            resp(i_tr).precue_event_onset_t_est*1000;
    end

    % main events
    % check for hits {'hit','miss','CR','FA_proper','FA'}
    resp(i_tr).button_presses_type = {}; %{'hit','miss','CR','FA_proper','FA'}
    resp(i_tr).button_presses_RT = []; % reaction time in ms (after event or closest to other event)
    resp(i_tr).event_response_type = repmat({'noevent'},1,max(p.stim.eventnum)); %{'hit','miss','CR','FA_proper','noevent'}
    resp(i_tr).event_response_RT = nan(1,max(p.stim.eventnum)); %reaction time or nan
    
    % first define response windows
    if any(~isnan(resp(i_tr).eventtype))
        t.respwin = (resp(i_tr).event_onset_times(~isnan(resp(i_tr).eventtype))+(p.targ_respwin/1000))*1000;
    end
    % loop across button presses
    for i_press = 1:numel(t.presses)
        if any(~isnan(resp(i_tr).eventtype))
            t.idx = t.presses(i_press)>=t.respwin(:,1) & t.presses(i_press)<=t.respwin(:,2);
            % not in any response window? --> FA; no reaction time
            if sum(t.idx)== 0 & ~t.idx_pre(i_press)
                resp(i_tr).button_presses_type{i_press} = 'FA';
                % check time to next event
                t.RT = t.presses(i_press) - (resp(i_tr).event_onset_times*1000);
                if any(t.RT>0) % bias to preceding events
                    resp(i_tr).button_presses_RT(i_press) = max(t.RT(t.RT>0));
                else % if no preceding event calculate negative RT
                    resp(i_tr).button_presses_RT(i_press) = max(t.RT);
                end
                resp(i_tr).button_presses_RT(i_press) = resp(i_tr).button_presses_RT(i_press);
                
            elseif t.idx_pre(i_press)
                % check whether button press is a pre-cue response is
                resp(i_tr).button_presses_type{i_press} = 'pre_cue_response';
            else
                if resp(i_tr).eventtype(t.idx)== 1
                    resp(i_tr).button_presses_type{i_press} = 'hit';
                else
                    resp(i_tr).button_presses_type{i_press} = 'FA_proper';
                end
                resp(i_tr).event_response_type{t.idx} = resp(i_tr).button_presses_type{i_press};
                % reaction time?
                resp(i_tr).button_presses_RT(i_press) = t.presses(i_press)-(resp(i_tr).event_onset_times(t.idx)*1000);
                resp(i_tr).event_response_RT(t.idx) = resp(i_tr).button_presses_RT(i_press) ;
            end
           
        else
            if t.idx_pre(i_press)
                resp(i_tr).button_presses_type{i_press} = 'pre_cue_response';
            else
                % no events --> all presses are false alarms; no reaction time
                resp(i_tr).button_presses_type{i_press} = 'FA';
                resp(i_tr).button_presses_RT(i_press) = nan;
            end
        end
    end
    % check for correct responses or misses or correct rejections
    for i_ev = 1:resp(i_tr).eventnum
        if strcmp(resp(i_tr).event_response_type{i_ev},'noevent')
            if resp(i_tr).eventtype(i_ev)==1
                resp(i_tr).event_response_type{i_ev}='miss';
            else
                resp(i_tr).event_response_type{i_ev}='CR';
            end
            resp(i_tr).event_response_RT(i_ev) = nan;
        end
    end



      

    % old wait routine
    % crttime2 = GetSecs;
    % t.timetowait = ((p.ITI(1)/1000)+RDKin.trial.duration)-(crttime2 - inittime);
    % ttt=WaitSecs(t.timetowait);
    % crttime3 = GetSecs; % for troubleshooting
    
    
end

% stop response recording
KbQueueRelease(ps.RespDev);

% send start trigger
if flag_training~=1
    PRPX_sendRecTrigger('stop')
end



end