function [] = pres_instruction(p,ps,RDK,blocknum,randmat,con_flag)
%PRES_INSTRUCTION present instructions
%   p           = parameters
%   ps          = screen parameters

%% initial parameters
INST.flag =                 {'training';'experiment'};
INST.text{1} =              [...                % text for training
    sprintf('TRAINING - block %1.0f',blocknum)...
    '\n\n\nBitte fixieren Sie kontinuierlich das Fixationskreuz!'...
    '\n\n\nReagieren Sie sobald einer der Balken des Kreuzes kürzer wird'...
    '\n\nmit einem Druck auf die Leertaste.'...
    '\n\n\nSobald das Fixationskreuz die Farbe ändert, achten Sie auf die'...
    '\n\nangezeigten Punktewolken. Wenn sich ein Teil der beachteten'...
    '\n\nPunkte kurz kohärent in eine Richtung bewegt, drücken Sie'...
    '\n\ndie Leertaste. Seien Sie dabei so schnell und akurat wie möglich.'...
    sprintf('\n\n\nNutzen Sie für die Antworten die %s Hand.', randmat.mats.responsehand{find(randmat.mats.block==blocknum,1,'first')})
    ];
INST.text{2} =              [...                % text for experiment
    sprintf('EXPERIMENT - Block %1.0f von %1.0f',blocknum,p.stim.blocknum)...
    '\n\n\nBitte fixieren Sie kontinuierlich das Fixationskreuz!'...
    '\n\n\nReagieren Sie sobald einer der Balken des Kreuzes kürzer wird'...
    '\n\nmit einem Druck auf die Leertaste.'...
    '\n\n\nSobald das Fixationskreuz die Farbe ändert, achten Sie auf die'...
    '\n\nangezeigten Punktewolken. Wenn sich ein Teil der beachteten'...
    '\n\nPunkte kurz kohärent in eine Richtung bewegt, drücken Sie'...
    '\n\ndie Leertaste. Seien Sie dabei so schnell und akurat wie möglich.'...
    sprintf('\n\n\nNutzen Sie für die Antworten die %s Hand.', randmat.mats.responsehand{find(randmat.mats.block==blocknum,1,'first')})
    ];
INST.text

% position of fixation cross
INST.crs.xCoords =          [-p.crs_dims(1) p.crs_dims(1) 0 0];
INST.crs.yCoords =          [0 0 -p.crs_dims(2) p.crs_dims(2)];
INST.crs.allCoords =        [INST.crs.xCoords; INST.crs.yCoords];
INST.crs.offset =           [0 180; 0 280]; % [x y] offset of cue cross

%% present text and dots
WaitSecs(0.5);
key.rkey=KbName('Space');                               % identify space bar
count_frame = 1;                                        % counter for presnetation of dot
fprintf(1,'\n%s - Instruktion',INST.flag{con_flag})
[key.keyisdown,key.secs,key.keycode] = KbCheck;         % initial check for preallocation

while ~(key.keycode(key.rkey)==1)                       % continuously present instruction
    % for i_fr = 1:2000
    [key.keyisdown,key.secs,key.keycode] = KbCheck;
    Screen('TextSize', ps.window, 15);
    DrawFormattedText(ps.window, INST.text{con_flag}, 'center', ps.screenYpixels * 0.1, p.stim_color);

    if con_flag ~=3
        % draw fixation cross   
        Screen('DrawLines', ps.window, INST.crs.allCoords, p.crs_width, p.cue_color(1,:), [ps.xCenter-INST.crs.offset(1,1) ps.yCenter+INST.crs.offset(1,2)], 2);
        DrawFormattedText(ps.window, 'Beachte links', 'center',  ps.yCenter+INST.crs.offset(1,2)+30, p.cue_color(1,:),[],[],[],[],[]);
        Screen('DrawLines', ps.window, INST.crs.allCoords, p.crs_width, p.cue_color(2,:), [ps.xCenter+INST.crs.offset(2,1) ps.yCenter+INST.crs.offset(2,2)], 2);
        DrawFormattedText(ps.window, 'Beachte rechts', 'center', ps.yCenter+INST.crs.offset(2,2)+30, p.cue_color(2,:),[],[],[],[],[]);
        
        % draw stimulus
        ovalOutRect = [0 0 p.stim_out_diam p.stim_out_diam];
        ovalInRect = [0 0 p.stim_in_diam p.stim_in_diam];
        maxDiameter = max(ovalOutRect) * 1.01;
        % shift rectangles to points and create matrix ([left outer; left inner; right outer; right inner])
        OvalsPos = [CenterRectOnPointd(ovalOutRect,ps.xCenter-p.stim_pos(1), ps.yCenter+p.stim_pos(2)+mean(INST.crs.offset(:,2)))'...
            CenterRectOnPointd(ovalInRect,ps.xCenter-p.stim_pos(1), ps.yCenter+p.stim_pos(2)+mean(INST.crs.offset(:,2)))'...
            CenterRectOnPointd(ovalOutRect,ps.xCenter+p.stim_pos(1), ps.yCenter+p.stim_pos(2)+mean(INST.crs.offset(:,2)))'...
            CenterRectOnPointd(ovalInRect,ps.xCenter+p.stim_pos(1), ps.yCenter+p.stim_pos(2)+mean(INST.crs.offset(:,2)))'];
        OvalCol = repmat(p.scr_color',1,4); OvalCol(:,1)= p.stim_color'; OvalCol(:,3)= p.stim_color';
        Screen('FillOval', ps.window, OvalCol, OvalsPos(:,[1 2 3 4]), maxDiameter);
        
        % draw target
        Screen('FillArc', ps.window, mean(p.targ_color,1), OvalsPos(:,1), 161, p.targ_gap(1));
        
        % draw inner circle
        Screen('FillOval', ps.window, OvalCol(:,[2 4]), OvalsPos(:,[2 4]), maxDiameter);
        
    end

    Screen('Flip', ps.window, 0);                   % flip screen
end

%% alert before experiment starts
Screen('TextSize', ps.window, 15);
DrawFormattedText(ps.window, 'Gleich geht es los.','center','center', p.stim_color);
Screen('Flip', ps.window, 0);
t=WaitSecs(2);


end

