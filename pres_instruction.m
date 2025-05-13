function [] = pres_instruction(p,ps,RDK,blocknum,randmat,con_flag, key)
%PRES_INSTRUCTION present instructions
%   p           = parameters
%   ps          = screen parameters

%% initial parameters
INST.flag =                 {'training';'experiment'};

% with pre-cue task
% INST.text{1} =              [...                % text for training
%     sprintf('TRAINING - block %1.0f',blocknum)...
%     '\n\n\nBitte fixieren Sie kontinuierlich das Fixationskreuz!'...
%     '\n\nReagieren Sie sobald einer der Balken des Kreuzes k�rzer wird'...
%     '\nmit einem Druck auf die Leertaste.'...
%     '\n\nSobald das Fixationskreuz die Farbe �ndert, achten Sie auf die'...
%     '\nangezeigten Punktewolken. Wenn sich ein Teil der beachteten'...
%     '\nPunkte kurz koh�rent in eine Richtung bewegt, dr�cken Sie'...
%     '\ndie Leertaste. Seien Sie dabei so schnell wie m�glich.'...
%     sprintf('\n\nNutzen Sie f�r die Antworten die %s Hand.', randmat.mats.responsehand{find(randmat.mats.block==1,1,'first')})
%     ];
% INST.text{2} =              [...                % text for experiment
%     sprintf('EXPERIMENT - Block %1.0f von %1.0f',blocknum,p.stim.blocknum)...
%     '\n\n\nBitte fixieren Sie kontinuierlich das Fixationskreuz!'...
%     '\n\nReagieren Sie sobald einer der Balken des Kreuzes k�rzer wird'...
%     '\nmit einem Druck auf die Leertaste.'...
%     '\n\nSobald das Fixationskreuz die Farbe �ndert, achten Sie auf die'...
%     '\nangezeigten Punktewolken. Wenn sich ein Teil der beachteten'...
%     '\nPunkte kurz koh�rent in eine Richtung bewegt, dr�cken Sie'...
%     '\ndie Leertaste. Seien Sie dabei so schnell wie m�glich.'...
%     sprintf('\n\nNutzen Sie f�r die Antworten die %s Hand.', randmat.mats.responsehand{find(randmat.mats.block==blocknum,1,'first')})
%     ];

% without pre-cue task
INST.text{1} =              [...                % text for training
    sprintf('TRAINING - block %1.0f',blocknum)...
    '\n\n\nBitte fixieren Sie kontinuierlich das Fixationskreuz!'...
    '\n\nSobald das Fixationskreuz die Farbe �ndert, achten Sie auf die'...
    '\nPunktewolken der angezeigten Farbe. Wenn sich ein Teil der beachteten'...
    '\nPunkte kurz koh�rent in eine Richtung bewegt, dr�cken Sie die Leertaste.'...
    '\nSeien Sie dabei so schnell wie m�glich.'...
    sprintf('\n\nNutzen Sie f�r die Antworten die %s Hand.', randmat.mats.responsehand{find(randmat.mats.block==1,1,'first')})
    ];
INST.text{2} =              [...                % text for experiment
    sprintf('EXPERIMENT - Block %1.0f von %1.0f',blocknum,p.stim.blocknum)...
    '\n\n\nBitte fixieren Sie kontinuierlich das Fixationskreuz!'...
    '\n\nSobald das Fixationskreuz die Farbe �ndert, achten Sie auf die'...
    '\nPunktewolken der angezeigten Farbe. Wenn sich ein Teil der beachteten'...
    '\nPunkte kurz koh�rent in eine Richtung bewegt, dr�cken Sie die Leertaste.'...
    '\nSeien Sie dabei so schnell wie m�glich.'...
    sprintf('\n\nNutzen Sie f�r die Antworten die %s Hand.', randmat.mats.responsehand{find(randmat.mats.block==blocknum,1,'first')})
    ];


% position of fixation cross
% INST.crs.xCoords =          [-p.crs_dims(1) p.crs_dims(1) 0 0];
% INST.crs.yCoords =          [0 0 -p.crs_dims(2) p.crs_dims(2)];
% INST.crs.allCoords =        [INST.crs.xCoords; INST.crs.yCoords];
% INST.crs.offset =           [0 180; 0 280]; % [x y] offset of cue cross

% pixels for shift into 4 quadrants
quadshift = [p.scr_res(1)*(1/4) p.scr_res(2)*(1/4); p.scr_res(1)*(3/4) p.scr_res(2)*(1/4); ...
    p.scr_res(1)*(1/4) p.scr_res(2)*(3/4); p.scr_res(1)*(3/4) p.scr_res(2)*(3/4)];

%% create fixation cross
ps.center = [ps.xCenter ps.yCenter];
p.crs.half = p.crs.dims/2;
p.crs.bars = [-p.crs.half p.crs.half 0 0; 0 0 -p.crs.half p.crs.half];

posshift1 = [repmat(-610,1,4); repmat(-200,1,4)];
posshift2 = [repmat(-610,1,4); repmat(-180,1,4)];
posshift3 = [repmat(-610,1,4); repmat(-160,1,4)];

p.crs.lines = {};
t.lines = [];
t.bars1 = p.crs.bars + posshift1 ;
t.bars2 = p.crs.bars + posshift2 ;
t.bars3 = p.crs.bars + posshift3 ;


%% present text in quadrants
% draw text and stimuli (before shifting to the quadrants)
% offscreen window
[ps.offwin,ps.offrect]=Screen('OpenOffscreenWindow',p.scr_num, [0 0 0 0], [0 0 p.scr_res(1)/2 p.scr_res(2)/2], [], [], []);
% get center of offscreen window
[ps.xCenter_off, ps.yCenter_off] = RectCenter(ps.offrect);

% draw instruction text
Screen('TextSize', ps.offwin, 14);
% DrawFormattedText(tx.instruct, INST.text{1}, tx.xCenter_off, p.scr_res(1)/2 * 0.1, p.stim_color);
DrawFormattedText(ps.offwin, INST.text{con_flag}, 'center', p.scr_res(1)/2 * 0.1, p.crs.color);

% draw text, how to start trial
DrawFormattedText(ps.offwin, 'Mit LEERTASTE geht es los!', 'center', p.scr_res(2)/2 * 0.85, p.crs.color);


% draw fixation cross with relevant cueing
Screen('DrawLines', ps.offwin, t.bars1, p.crs.width, RDK.RDK(1).col(1,:), ps.center, 0);
Screen('DrawLines', ps.offwin, t.bars2, p.crs.width, RDK.RDK(2).col(1,:), ps.center, 0);
Screen('DrawLines', ps.offwin, t.bars3, p.crs.width, RDK.RDK(3).col(1,:), ps.center, 0);

% draw text for instruction (which RDK to attend)
DrawFormattedText(ps.offwin, 'RDKs zu beachten:', p.scr_res(1)/2 * 0.41, p.scr_res(1)/2 * 0.359, p.crs.color);
DrawFormattedText(ps.offwin, 'RDKs zu beachten:', p.scr_res(1)/2 * 0.41, p.scr_res(1)/2 * 0.381, p.crs.color);
DrawFormattedText(ps.offwin, 'RDKs zu beachten:', p.scr_res(1)/2 * 0.41, p.scr_res(1)/2 * 0.402, p.crs.color);

DrawFormattedText(ps.offwin, RDK.RDK(1).col_label, p.scr_res(1)/2 * 0.57, p.scr_res(1)/2 * 0.359, RDK.RDK(1).col(1,:));
DrawFormattedText(ps.offwin, RDK.RDK(2).col_label, p.scr_res(1)/2 * 0.57, p.scr_res(1)/2 * 0.381, RDK.RDK(2).col(1,:));
DrawFormattedText(ps.offwin, RDK.RDK(3).col_label, p.scr_res(1)/2 * 0.57, p.scr_res(1)/2 * 0.402, RDK.RDK(3).col(1,:));

DrawFormattedText(ps.offwin, RDK.RDK(2).col_label, p.scr_res(1)/2 * 0.61, p.scr_res(1)/2 * 0.359, RDK.RDK(2).col(1,:));
DrawFormattedText(ps.offwin, RDK.RDK(3).col_label, p.scr_res(1)/2 * 0.61, p.scr_res(1)/2 * 0.381, RDK.RDK(3).col(1,:));
DrawFormattedText(ps.offwin, RDK.RDK(1).col_label, p.scr_res(1)/2 * 0.61, p.scr_res(1)/2 * 0.402, RDK.RDK(1).col(1,:));


for i_quad = 1:4 % shifst to quadrants
    newpos_stim(:,i_quad) = ...
        CenterRectOnPointd(ps.offrect,quadshift(i_quad,1),quadshift(i_quad,2))';
end

% show info outside
fprintf('\nVersuchsperson startet %s Block %1.0f mit Leertaste...', INST.flag{con_flag},blocknum)

% [key.pressed, key.firstPress]=KbQueueCheck;
key.rkey=key.SPACE;
[key.keyisdown,key.secs,key.keycode] = KbCheck; 
while ~(key.keycode(key.rkey)==1)                       % continuously present feedback (wait for q)
    [key.keyisdown,key.secs,key.keycode] = KbCheck;
    % draw text
    Screen('DrawTextures', ps.window, repmat(ps.offwin,1,4),[], newpos_stim, [], [], [], []);
    Screen('Flip', ps.window, 0);                   % flip screen
end
try Screen('Close', ps.offwin); end



%% alert before experiment starts
[ps.offwin,ps.offrect]=Screen('OpenOffscreenWindow',p.scr_num, [0 0 0 0], [0 0 p.scr_res(1)/2 p.scr_res(2)/2], [], [], []);
% get center of offscreen window
[ps.xCenter_off, ps.yCenter_off] = RectCenter(ps.offrect);

% draw instruction text
Screen('TextSize', ps.offwin, 14);
DrawFormattedText(ps.offwin, 'Gleich geht es los.','center','center', p.crs.color);
for i_quad = 1:4 % shifst to quadrants
    newpos_stim(:,i_quad) = ...
        CenterRectOnPointd(ps.offrect,quadshift(i_quad,1),quadshift(i_quad,2))';
end
Screen('DrawTextures', ps.window, repmat(ps.offwin,1,4),[], newpos_stim, [], [], [], []);
Screen('Flip', ps.window, 0);
fprintf('Block beginnt!\n')
t=WaitSecs(2);
try Screen('Close', ps.offwin); end

end

