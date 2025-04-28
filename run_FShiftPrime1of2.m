function [] = run_FShiftBase(sub,flag_training, flag_isolum, flag_block)
% run_FShiftBase(sub,flag_training, flag_isolum, flag_block)
%   runs experiment SSVEP_FShiftBase
%       sub:            participant number
%       flag_training:  1 = do training
%       flag_isolum:    1 = do isoluminance adjustment
%       flag_block:     1 = start with block 1
%           e.g. run_FShiftBase(1,1, 0, 1)



% Maria Dotzer, Christopher Gundlach, Leipzig, 2020

if nargin < 4
    help run_FShiftBase
    return
end

%% parameters
% sub = 1; flag_training = 1; flag_block = 1; flag_isolum = 1;
% design
p.sub                   = sub;                  % subject number
p.flag_block            = flag_block;           % block number to start
p.flag_training         = flag_training;        % do training

p.ITI                   = [1000 1000];          % inter trial interval in ms
p.targ_respwin          = [200 1000];           % time window for responses in ms

% screen
p.scr_num               = 1;                    % screen number
p.scr_res               = [1920 1080];          % resolution
p.scr_refrate           = 480;                  % refresh rate in Hz (e.g. 85)
p.scr_color             = [0.05 0.05 0.05 1];      % default: [0.05 0.05 0.05 1]; ; color of screen [R G B Alpha]
p.scr_imgmultipl        = 4;

% some isoluminace parameters
p.isol.TrlAdj           = 5;                    % number of trials used for isoluminance adjustment
p.isol.MaxStd           = 10;                   % standard deviation tolerated
p.isol.run              = false;                % isoluminance run?
p.isol.override         = [];                   % manually set colors for RDK1 to RDKXs e.g. []
% p.isol.override         = [0 0.0627450980392157 0.156862745098039 1;0 0.0823529411764706 0 1;0.0862745098039216 0.0345098039215686 0 1];
% p.isol.override         = [0 0.332549019607843 0.831372549019608 1;0 0.439215686274510 0 1;0.454901960784314 0.181960784313726 0 1];
p.isol.bckgr            = p.scr_color(1:3)+0.2;          % isoluminant to background or different color?
% p.isol.bckgr            = p.scr_color;          % isoluminant to background or different color?


% stimplan
p.stim.condition        = [1 2];              % conditions
p.stim.RDK2attend       = [1 2];            % defines which RDK to attend in which condition
p.stim.eventnum         = [0 0 1 2];        % ratio of eventnumbers
p.stim.con_repeats_e    = 60;               % trial number/repeats for each eventnum in experiment
p.stim.con_repeats_t    = 2;                % trial number/repeats for each eventnum in training
p.stim.time_postcue     = 2;                % post.cue time in s
p.stim.time_precue      = [1.5 2];          % precue time in s; [upper lower] for randomization
p.stim.event.type       = 2;                % types of events (1 = targets only, 2 = targets + distrators)
p.stim.event.length     = 0.3;              % lengt of events in s
p.stim.event.min_onset  = 0.2;              % min post-cue time before event onset in s
p.stim.event.min_offset = 0;                % min offset from target end to end of trial in s
p.stim.event.min_dist   = 0.8;              % min time between events in s
p.stim.blocknum         = 8;                % number of blocks
p.stim.ITI              = [1 1];            % ITI range in seconds
p.stim.frames_postcue   = p.stim.time_postcue*p.scr_refrate;

% stimuli
RDK.RDK(1).size             = [360 360];                % width and height of RDK in pixel; only even values 
RDK.RDK(1).centershift      = [0 0];                    % position of RDK center; x and y deviation from center in pixel
RDK.RDK(1).col              = [0 1 0 1; p.scr_color(1:3) 0];% "on" and "off" color
RDK.RDK(1).freq             = 30;                       % flicker frequency, frequency of a full "on"-"off"-cycle
RDK.RDK(1).mov_freq         = 120;                      % Defines how frequently the dot position is updated; 0 will adjust the update-frequency to your flicker frequency (i.e. dot position will be updated with every "on"-and every "off"-frame); 120 will update the position for every frame for 120Hz or for every 1. quadrant for 480Hz 
RDK.RDK(1).num              = 85;                      % number of dots
RDK.RDK(1).mov_speed        = 1;                        % movement speed in pixel
RDK.RDK(1).mov_dir          = [0 1; 0 -1; -1 0; 1 0];   % movement direction  [0 1; 0 -1; -1 0; 1 0] = up, down, left, right
RDK.RDK(1).dot_size         = 12;
 
RDK.RDK(2).size             = RDK.RDK(1).size;                
RDK.RDK(2).centershift      = RDK.RDK(1).centershift;                  
RDK.RDK(2).col              = [0 0.4 1 1; p.scr_color(1:3) 0];
RDK.RDK(2).freq             = 24;
RDK.RDK(2).mov_freq         = RDK.RDK(1).mov_freq;
RDK.RDK(2).num              = RDK.RDK(1).num;
RDK.RDK(2).mov_speed        = RDK.RDK(1).mov_speed;
RDK.RDK(2).mov_dir          = RDK.RDK(1).mov_dir;
RDK.RDK(2).dot_size         = RDK.RDK(1).dot_size;

RDK.RDK(3).size             = RDK.RDK(1).size;              
RDK.RDK(3).centershift      = RDK.RDK(1).centershift;                   
RDK.RDK(3).col              = [1 0.4 0 1; p.scr_color(1:3) 0];
RDK.RDK(3).freq             = 20;
RDK.RDK(3).mov_freq         = RDK.RDK(1).mov_freq;
RDK.RDK(3).num              = RDK.RDK(1).num;
RDK.RDK(3).mov_speed        = RDK.RDK(1).mov_speed;
RDK.RDK(3).mov_dir          = RDK.RDK(1).mov_dir;
RDK.RDK(3).dot_size         = RDK.RDK(1).dot_size;

RDK.event.type              = 'globalmotion';               % event type global motion
RDK.event.duration          = p.stim.event.length;          % time of coherent motion
RDK.event.coherence         = .4;                            % percentage of coherently moving dots 0.7
RDK.event.direction         = RDK.RDK(1).mov_dir;           % movement directions for events

% fixation cross
p.crs.color                 = [0.8 0.8 0.8 1];    % color of fixation cross
p.crs.dims                  = [16];             % dimension of fixation cross
p.crs.width                 = 2;                % width of fixation cross
p.crs.cutout                = 0;                % 1 = no dots close to fixation cross

% trigger
p.trig.rec_start        = 253;                  % trigger to start recording
p.trig.rec_stop         = 254;                  % trigger to stop recording
p.trig.tr_start         = 11;                   % trial start; main experiment
p.trig.tr_stop          = 13;                   % trial end; main experiment
p.trig.tr_cue_type      = [100 200];            % cue type: [RDK1 RDK2]
p.trig.type             = [10 20; 1 2];         % [first: target, distractor; second: target, distractor]
p.trig.button           = 60;                   % button press
p.trig.event_type       = [80 90];              % target, distractor
p.trig.event_dir        = 1:size(RDK.RDK(1).mov_dir); % movement direction


% logfiles
p.log.path              = '/home/pc/matlab/user/christopher/SSVEP_FShiftBase/logfiles/';
p.log.exp_name          = 'SSVEP_FShiftBase';
p.log.add               = '_a';


%% check for logfile being present
filecheck=dir(sprintf('%sVP%02.0f_timing*',p.log.path,p.sub));
if ~isempty(filecheck)
    reply = input(sprintf('\nVP%02.0f existiert bereits. Datei überschreiben? [j/n]... ',p.sub),'s');
    if strcmp(reply,'j')
        p.filename = sprintf('VP%02.0f_timing',p.sub);
    else
        [temp name_ind]=max(cellfun(@(x) numel(x), {filecheck.name}));
        p.filename = sprintf('%s%s',filecheck(name_ind).name(1:end-4),p.log.add);
    end
else
    p.filename = sprintf('VP%02.0f_timing',p.sub);
end

t.isol = {};
% routine to check for older isoluminance adjustments
for i_file = 1:numel(filecheck)
    t.in = load(fullfile(filecheck(i_file).folder,filecheck(i_file).name));
    t.datenum{i_file} = filecheck(i_file).datenum;
    t.isol{i_file} = t.in.p.isol;
    
end



%% Screen init
ps.input = struct('ScrNum',p.scr_num,'RefRate',p.scr_refrate,'PRPXres',p.scr_res,'BckGrCol',p.scr_color,'PRPXmode',2);
[~, ps.screensize, ps.xCenter, ps.yCenter, ps.window, ps.framerate, ps.RespDev, ps.keymap] = PTExpInit_GLSL(ps.input,1);

% some initial calculations
% fixation cross
ps.center = [ps.xCenter ps.yCenter];
p.crs.half = p.crs.dims/2;
p.crs.bars = [-p.crs.half p.crs.half 0 0; 0 0 -p.crs.half p.crs.half];

% shift into 4 quadrants (running with 480 Hz)
ps.shift = [-ps.xCenter/2, -ps.yCenter/2; ps.xCenter/2, -ps.yCenter/2;... % shifts to four quadrants: upper left, upper right, lower left, lower right
    -ps.xCenter/2, ps.yCenter/2; ps.xCenter/2, ps.yCenter/2];

p.crs.lines = [];
for i_quad=1:p.scr_imgmultipl
    p.crs.lines = cat(2, p.crs.lines, [p.crs.bars(1,:)+ps.shift(i_quad,1); p.crs.bars(2,:)+ps.shift(i_quad,2)]); %array with start and end points for the fixation cross lines, for all four quadrants
end

%% keyboard and ports setup ???
KbName('UnifyKeyNames')
Buttons = [KbName('ESCAPE') KbName('Q') KbName('SPACE') KbName('j') KbName('n') KbName('1!') KbName('2@') KbName('3#')];
RestrictKeysForKbCheck(Buttons);
key.keymap=false(1,256);
key.keymap(Buttons) = true;
key.keymap_ind = find(key.keymap);
[key.ESC, key.SECRET, key.SPACE, key.YES, key.NO] = deal(...
    Buttons(1),Buttons(2),Buttons(3),Buttons(4),Buttons(5));

%% start experiment
% initialize randomization of stimulation frequencies and RDK colors [RDK1 and RDK2 task relevant RDK3 not]
rand('state',p.sub)
[RDK.RDK(:).col_init] = deal(RDK.RDK(:).col);
[RDK.RDK(:).freq_init] = deal(RDK.RDK(:).freq);
[RDK.RDK(:).col] = deal(RDK.RDK(randperm(3)).col);
[RDK.RDK(:).freq] = deal(RDK.RDK(randperm(3)).freq);

% initialize blank variables
timing = []; button_presses = []; resp = []; randmat = [];

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
        rand('state',p.sub*i_bl) % determine randstate
        randmat.training{i_bl} = rand_FShiftBase(p, RDK,  1);
        [timing.training{i_bl},button_presses.training{i_bl},resp.training{i_bl}] = ...
            pres_FShiftBase(p, ps, key, RDK, randmat.training{i_bl}, i_bl,1);
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

%% then isoluminance adjustment
% do the heterochromatic flicker photometry
if flag_isolum == 1
%     
%     PsychDefaultSetup(2);
%     Datapixx('Open');
%     Datapixx('SetPropixxDlpSequenceProgram', 0);
%     Datapixx('RegWrRd');
     
    % initial colors
    p.isol.init_cols = cell2mat({RDK.RDK.col}'); p.isol.init_cols = p.isol.init_cols(1:2:end,1:3);
    
    % start isoluminance script only RGB output (no alpha)
    [Col2Use] = PRPX_IsolCol_480_adj([p.isol.bckgr(1:3); p.isol.init_cols],p.isol.TrlAdj,p.isol.MaxStd,0,RDK.RDK(1).size);
%     [Col2Use] = PRPX_IsolCol_adj([p.isol.bckgr(1:3); p.isol.init_cols],p.isol.TrlAdj,p.isol.MaxStd,0,RDK.RDK(1).size);
    
    for i_RDK = 1:numel(RDK.RDK)
        RDK.RDK(i_RDK).col(1,:) = [Col2Use(1+i_RDK,:) 1];
    end
    % index function execution
    p.isol.run = sprintf('originally run: %s',datestr(now));
    p.isol.coladj = [Col2Use(2:end,:) ones(size(Col2Use,1)-1,1)];
    save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
    
    fprintf('\nadjusted colors:\n')
    for i_col = 1:size(p.isol.coladj,1)
        fprintf('RDK%1.0f [%1.4f %1.4f %1.4f %1.4f]\n', i_col,p.isol.coladj(i_col,:))
    end
    
    Screen('CloseAll')
    Datapixx('SetPropixxDlpSequenceProgram', 0);
    Datapixx('RegWrRd');
    Datapixx('close');
else
    % select colors differently
    fprintf(1,'\nKeine Isoluminanzeinstellung. Wie soll verfahren werden?')
    % specify options
    % option1: use default values
    isol.opt(1).available = true;
    t.cols = cell2mat({RDK.RDK(:).col}');
    isol.opt(1).colors = t.cols(1:2:end,:);
    isol.opt(1).text = sprintf('default: [%1.2f %1.2f %1.2f %1.2f] [%1.2f %1.2f %1.2f %1.2f] [%1.2f %1.2f %1.2f %1.2f]',isol.opt(1).colors');
    % option2: use isoluminance values of previously saved dataset
    if ~isempty(t.isol) % file loaded?
        [t.t t.idx] = max(cell2mat(t.datenum));
        isol.opt(2).available = true;
        isol.opt(2).colors = t.isol{t.idx}.coladj(1:end,:);
        isol.opt(2).text = sprintf('aus gespeicherter Datei: [%1.2f %1.2f %1.2f %1.2f] [%1.2f %1.2f %1.2f %1.2f] [%1.2f %1.2f %1.2f %1.2f]',isol.opt(2).colors');
    else
        isol.opt(2).available = false;
        isol.opt(2).colors = [];
        isol.opt(2).text = [];
    end
    % option3: use manual override
    if ~isempty(p.isol.override)
        isol.opt(3).available = true;
        isol.opt(3).colors = p.isol.override;
        isol.opt(3).text = sprintf('manuell definiert in p.isol override: [%1.2f %1.2f %1.2f %1.2f] [%1.2f %1.2f %1.2f %1.2f] [%1.2f %1.2f %1.2f %1.2f]',isol.opt(3).colors');
    else
        isol.opt(3).available = false;
        isol.opt(3).colors = [];
        isol.opt(3).text = [];
    end
    % check for buttons
    IsoButtons = Buttons(6:8);
    isol.prompt.idx = find([isol.opt(:).available]);
    t.prompt = [];
    for i_prompt = 1:numel(isol.prompt.idx)
        t.prompt = [t.prompt sprintf('\n(%1.0f) %s',i_prompt,isol.opt(isol.prompt.idx(i_prompt)).text)];
    end
    
    % display options
    fprintf('%s',t.prompt)
    inp.prompt_check = 0;
    while inp.prompt_check == 0             % loop to check for correct input
        [key.keyisdown,key.secs,key.keycode] = KbCheck;
        if any(key.keycode)
            inp.prompt_check = 1;
        end
        Screen('Flip', ps.window, 0);
    end
    Col2Use = isol.opt(isol.prompt.idx(key.keycode(IsoButtons(1:numel(isol.prompt.idx)))==1)).colors;
    % use selected colors
    for i_RDK = 1:numel(RDK.RDK)
        RDK.RDK(i_RDK).col(1,:) = Col2Use(i_RDK,:);
    end
    % index function execution
    switch isol.prompt.idx(key.keycode(IsoButtons(1:numel(isol.prompt.idx)))==1)
        case 1
            p.isol.run = sprintf('default at %s',datestr(now));
        case 2
            p.isol.run = sprintf('reloaded at %s from %s',datestr(now),datestr(t.datenum{t.idx}));
        case 3
            p.isol.run = sprintf('override at %s',datestr(now));
    end
    p.isol.coladj = Col2Use;
%     save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
    
    fprintf('\nselected colors:\n')
    for i_col = 1:size(p.isol.coladj,1)
        fprintf('RDK%1.0f [%1.4f %1.4f %1.4f]\n', i_col,p.isol.coladj(i_col,:))
    end
end

%% redo initialization
ps.input = struct('ScrNum',p.scr_num,'RefRate',p.scr_refrate,'PRPXres',p.scr_res,'BckGrCol',p.scr_color,'PRPXmode',2);
[~, ps.screensize, ps.xCenter, ps.yCenter, ps.window, ps.framerate, ps.RespDev, ps.keymap] = PTExpInit_GLSL(ps.input,1);

% some initial calculations
% fixation cross
ps.center = [ps.xCenter ps.yCenter];
p.crs.half = p.crs.dims/2;
p.crs.bars = [-p.crs.half p.crs.half 0 0; 0 0 -p.crs.half p.crs.half];

% shift into 4 quadrants (running with 480 Hz)
ps.shift = [-ps.xCenter/2, -ps.yCenter/2; ps.xCenter/2, -ps.yCenter/2;... % shifts to four quadrants: upper left, upper right, lower left, lower right
    -ps.xCenter/2, ps.yCenter/2; ps.xCenter/2, ps.yCenter/2];

p.crs.lines = [];
for i_quad=1:p.scr_imgmultipl
    p.crs.lines = cat(2, p.crs.lines, [p.crs.bars(1,:)+ps.shift(i_quad,1); p.crs.bars(2,:)+ps.shift(i_quad,2)]); %array with start and end points for the fixation cross lines, for all four quadrants
end

% keyboard setup
KbName('UnifyKeyNames')
Buttons = [KbName('ESCAPE') KbName('Q') KbName('SPACE') KbName('j') KbName('n') KbName('1!') KbName('2@') KbName('3#')];
RestrictKeysForKbCheck(Buttons);
key.keymap=false(1,256);
key.keymap(Buttons) = true;
key.keymap_ind = find(key.keymap);
[key.ESC, key.SECRET, key.SPACE, key.YES, key.NO] = deal(...
    Buttons(1),Buttons(2),Buttons(3),Buttons(4),Buttons(5));


%% do training again?
% loop for training to be repeated
fprintf(1,'\nTraing starten (j/n)')
inp.prompt_check = 0;
while inp.prompt_check == 0             % loop to check for correct input
    [key.keyisdown,key.secs,key.keycode] = KbCheck;
    if key.keycode(key.YES)==1
        flag_trainend = 0; inp.prompt_check = 1;
    elseif key.keycode(key.NO)==1
        flag_trainend = 1; inp.prompt_check = 1;
    end
    Screen('Flip', ps.window, 0);
end

if ~exist('i_bl'); i_bl = 1; end
while flag_trainend == 0 % do training until ended
    rand('state',p.sub*i_bl) % determine randstate
    randmat.training{i_bl} = rand_FShiftBase(p, RDK,  1);
    [timing.training{i_bl},button_presses.training{i_bl},resp.training{i_bl}] = ...
        pres_FShiftBase(p, ps, key, RDK, randmat.training{i_bl}, i_bl,1);
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


%% present each block
% randomization
rand('state',p.sub);                         % determine randstate
randmat.experiment = rand_FShiftBase(p, RDK,  0);    % randomization
for i_bl = p.flag_block:p.stim.blocknum
    % start experiment
    [timing.experiment{i_bl},button_presses.experiment{i_bl},resp.experiment{i_bl}] = ...
        pres_FShiftBase(p, ps, key, RDK, randmat.experiment, i_bl,0);
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

