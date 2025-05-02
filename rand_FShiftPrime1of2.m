function [ conmat ] = rand_FShiftPrime1of2(p,RDK,flag_training)
%rand_FShiftBase randomizes experimental conditions
% move onset only works for constant frequency for all RDKs (i.e. 120)




% set trial number etc
if flag_training~=0
    conmat.totaltrials = numel(p.stim.condition)*numel(p.stim.eventnum)*p.stim.con_repeats_t;
    conmat.totalblocks = 1;
else
    conmat.totaltrials = numel(p.stim.condition)*numel(p.stim.eventnum)*p.stim.con_repeats_e;
    conmat.totalblocks = p.stim.blocknum;
end
conmat.trialsperblock = conmat.totaltrials/conmat.totalblocks;



%% there is a bug here....onsets and movements don*t match



% matrix with onset times of on framesfor RDKs
[t.onframesonset, t.onframesonset_times]= deal(nan(numel(RDK.RDK),p.scr_refrate*max(p.stim.time_postcue)));
for i_rdk = 1:numel(RDK.RDK)
    t.mat = ceil(1:p.scr_refrate/RDK.RDK(i_rdk).freq:size(t.onframesonset,2));
    t.onframesonset(i_rdk,t.mat)=1;
    t.onframesonset_times(i_rdk,t.mat)=t.mat./p.scr_refrate;
end
% figure; plot(~isnan(t.onframesonset_times)')

% move
[t.movonset_frames, t.movonset_times]=deal(nan(1,p.scr_refrate*max(p.stim.time_postcue)));
t.mat = 1:p.scr_refrate/RDK.RDK(1).mov_freq:size(t.movonset_frames,2);
t.movonset_frames(t.mat)=1;
t.movonset_times(t.mat)=t.mat./p.scr_refrate;

[t.onframesonset2, t.onframesonset_times2]= deal(nan(numel(RDK.RDK),p.scr_refrate*max(p.stim.time_postcue)));
for i_rdk = 1:numel(RDK.RDK)
    t.movonset_idx = find(~isnan(t.movonset_frames))';
    t.idx = dsearchn(find(~isnan(t.movonset_frames))',find(~isnan(t.onframesonset(i_rdk,:)))');
    t.onframesonset2(i_rdk,t.movonset_idx(t.idx))=1;
end
%% start randomization
% randomize cue (1 = attend to RDK 1,2; 2 = attend to RDK 2,3; 3 = attend to RDK 3,1)
conmat.mats.cue = reshape(repmat(p.stim.condition,conmat.totaltrials/numel(p.stim.condition),1),1,[]);

% which RDKs are attended?
conmat.mats.RDK2attend = p.stim.RDK2attend(conmat.mats.cue,:)';

% randomize event numbers per trial across condition
conmat.mats.eventnum = nan(size(conmat.mats.cue));
for i_cue = 1:numel(p.stim.condition)
    % index all trials of cue condition
    t.idx = conmat.mats.cue==p.stim.condition(i_cue);
    conmat.mats.eventnum(t.idx) = [ ... % append + randsample
        repmat(p.stim.eventnum,1,floor(sum(t.idx)/numel(p.stim.eventnum))), ...
        randsample(p.stim.eventnum, mod(sum(t.idx),numel(p.stim.eventnum)))];
end

% randomize eventtype (1 = target; 2 = distractor)
conmat.mats.eventtype = nan(max(p.stim.eventnum),conmat.totaltrials);
t.evnum = unique(p.stim.eventnum(p.stim.eventnum>0));
for i_cue = 1:numel(p.stim.condition)
    % index all trials of cue condition
    t.idx1 = conmat.mats.cue==p.stim.condition(i_cue);
    for i_evnum = 1:numel(t.evnum)
        % index all events of eventnum
        t.idx2 = conmat.mats.eventnum==t.evnum(i_evnum);
        % what type of events can be present as event 1,2,3,4,...
        t.evtype = repmat(cell2mat(arrayfun(@(x,y) repmat(x,1,y),[1,2],p.stim.event.ratio,'UniformOutput',false)), ...% define events acc. to specified ratio
            t.evnum(i_evnum),1);

        % find all combinations
        % Convert each row into a cell array
        t.rowCells = mat2cell(t.evtype, ones(1, t.evnum(i_evnum)), size(t.evtype,2));

        % Now use ndgrid to create a grid for all combinations
        t.grids = {};
        [t.grids{1:t.evnum(i_evnum)}] = ndgrid(t.rowCells{:});

        % Stack the results
        t.evtype_comb = cell2mat( cellfun(@(x) x(:), t.grids, 'UniformOutput', false) );
        t.evtype_comb = t.evtype_comb.'; % optional transpose to match combvec output

        % now write into conmat.mats.eventtype
        conmat.mats.eventtype(1:(t.evnum(i_evnum)),t.idx1&t.idx2) = [ ...
            repmat(t.evtype_comb,1,floor(sum(t.idx1&t.idx2)/size(t.evtype_comb,2))), ...
            t.evtype_comb(:,randsample(size(t.evtype_comb,2), mod(sum(t.idx1&t.idx2),size(t.evtype_comb,2))))];
    end
end
% some checking
% sum(conmat.mats.cue==1 & conmat.mats.eventtype(1,:)==1 & conmat.mats.eventtype(2,:) == 1)

% determine event RDK and distribute evenly
conmat.mats.eventRDK = nan(max(p.stim.eventnum),conmat.totaltrials);
t.evnum = unique(p.stim.eventnum(p.stim.eventnum>0));
t.eventtype = conmat.mats.eventtype; t.eventtype(isnan(t.eventtype))=0;
t.uniqueseqs = unique(t.eventtype','rows')';
% see where we find the unique t.uniqueseqs in conmat.mats.eventtype
t.uniqueseqs_exp = permute(t.uniqueseqs,[1 2 3]); % expand the 
t.eventtype_exp = permute(t.eventtype,[1 3 2]); % expand the
t.match = (t.eventtype_exp==t.uniqueseqs); % look for the tuple at each position
[t.rowmatch, t.colmatch]=find(squeeze(all(t.match,1))); % only columns where all rows match are kept; squeeze all of it
% loop across relevant variables
for i_cue = 1:numel(p.stim.condition)
    % index all trials of cue condition
    t.idx1 = conmat.mats.cue==p.stim.condition(i_cue);
    % index the relevant RDKs for the specific condition
    t.evRDK = {p.stim.RDK2attend(i_cue,:),find(~ismember(1:numel(RDK.RDK),p.stim.RDK2attend(i_cue,:)))};
    for i_evseqs = 1:size(t.uniqueseqs,2) % across all sequencess of eventtypes
        t.idx3 = t.rowmatch'==i_evseqs;
        % bring all together
        t.idxall = t.idx1&t.idx3;
        if any(t.idxall&conmat.mats.eventnum~=0) % any trial
            % which RDKs are targets and which distractors
            t.rdkidx = t.uniqueseqs(:,i_evseqs)~=0; % find the RDKs relevant for this cue
            t.rdkidx(~t.rdkidx)=[];
            t.rdkcombs = CombVec(t.evRDK{1,double([t.uniqueseqs(t.rdkidx,i_evseqs)'])}); % find combinations of RDKs (across potentially two numbers of events)
            if size(t.rdkcombs,1) == 1
                t.rdkcombs(2,:)=nan;
            end
            % distribute these across the conditions
            conmat.mats.eventRDK(:,t.idxall) = [ ... repmat and resample if necessary
                repmat(t.rdkcombs,1,floor(sum(t.idxall)/size(t.rdkcombs,2))), ...
                t.rdkcombs(:,randsample(size(t.rdkcombs,2),mod(sum(t.idxall),size(t.rdkcombs,2))))];
        end

    end
end



% randomize event directions (according to RDK.event.direction)
conmat.mats.eventdirection = nan(max(p.stim.eventnum),conmat.totaltrials);
conmat.mats.eventdirection(~isnan(conmat.mats.eventtype)) = ...
    randsample(size(RDK.event.direction,1),sum(~isnan(conmat.mats.eventtype),"all"),true);

% pre-allocate possible presentation times (works only for two events)
% sync it roughly with RDK onset frames
conmat.mats.event_onset_frames = nan(max(p.stim.eventnum),conmat.totaltrials);
t.poss_frames = repmat( ...
    p.stim.event.min_onset<t.movonset_times & ...
    t.movonset_times<(p.stim.time_postcue-p.stim.event.length-p.stim.event.min_offset), ...
    numel(RDK.RDK),1) & ~isnan(t.onframesonset2);
t.poss_frames_1 = repmat( ...
    p.stim.event.min_onset<t.movonset_times & ...
    t.movonset_times<(p.stim.time_postcue-p.stim.event.length-p.stim.event.min_offset-p.stim.event.min_dist-0.1), ...
    numel(RDK.RDK),1) & ~isnan(t.onframesonset2);
t.poss_frames_2 = repmat( ...
    p.stim.event.min_onset+p.stim.event.min_dist<t.movonset_times & ...
    t.movonset_times<(p.stim.time_postcue-p.stim.event.length-p.stim.event.min_offset), ...
    numel(RDK.RDK),1) & ~isnan(t.onframesonset2);

% figure; plot(t.poss_frames_1')
% figure; plot(t.onframesonset')

% now loop across conditions, event type and RDK (for first position)
for i_con = 1:numel(p.stim.condition)
    % index all trials of cue condition
    t.idx1 = conmat.mats.cue==p.stim.condition(i_con);
    % index RDKs
    t.evRDK = {p.stim.RDK2attend(i_con,:),find(~ismember(1:numel(RDK.RDK),p.stim.RDK2attend(i_con,:)))};
    for i_evtype = 1:2 % target vs distractor
        % index event types
        t.idx2 = conmat.mats.eventtype==i_evtype;
        t.evRDKi = t.evRDK{i_evtype};
        for i_rdk = 1:numel(t.evRDKi)
            t.idx3 = conmat.mats.eventRDK == t.evRDKi(i_rdk);
            
            % for single events first
            t.idx4 = conmat.mats.eventnum == 1;
            t.idxall = t.idx1 & t.idx2(1,:) & t.idx3(1,:) & t.idx4;
            % possible time frames for the respective rdk
            t.idxframes = find(t.poss_frames(t.evRDKi(i_rdk),:));
            % if there are more events than possible positions
            if sum(t.idxall)>numel(t.idxframes)
                t.poss_frames_mat = [ ...
                    repmat(t.idxframes,1,floor(sum(t.idxall)/numel(t.idxframes))), ...
                    t.idxframes(randsample(numel(t.idxframes),mod(sum(t.idxall),numel(t.idxframes))))];
            else
                % create some kind of bins to have events evenly distributed and the assign them accordingly
                t.poss_frames_mat = [];
                % first find greates common divisor:
                t.binnum = gcd(sum(t.idxall),length(t.idxframes));
                t.binLabels = discretize(1:length(t.idxframes), linspace(1, length(t.idxframes)+1, t.binnum+1));
                % now in a loop randomly sample from bins
                t.poss_frames_mat = cell2mat(arrayfun(@(x) ...
                    t.idxframes(randsample(find(t.binLabels==x),sum(t.idxall)/t.binnum)), ...
                    unique(t.binLabels), 'UniformOutput', false));
                conmat.mats.event_onset_frames(1,t.idxall) = t.poss_frames_mat(randperm(numel(t.poss_frames_mat)));
            end

            % for double events
            % randomize first position
            t.idx4 = conmat.mats.eventnum == 2;
            t.idxall = t.idx1 & t.idx2(1,:) & t.idx3(1,:) & t.idx4;
            % possible time frames for the respective rdk
            t.idxframes = find(t.poss_frames_1(t.evRDKi(i_rdk),:));
            % if there are more events than possible positions
            if sum(t.idxall)>numel(t.idxframes)
                t.poss_frames_mat = [ ...
                    repmat(t.idxframes,1,floor(sum(t.idxall)/numel(t.idxframes))), ...
                    t.idxframes(randsample(numel(t.idxframes),mod(sum(t.idxall),numel(t.idxframes))))];
            else
                % create some kind of bins to have events evenly distributed and the assign them accordingly
                t.poss_frames_mat = [];
                % first find greatest common divisor:
                t.binnum = gcd(sum(t.idxall),length(t.idxframes));
                t.binLabels = discretize(1:length(t.idxframes), linspace(1, length(t.idxframes)+1, t.binnum+1));
                % now in a loop randomly sample from bins
                t.poss_frames_mat = cell2mat(arrayfun(@(x) ...
                    t.idxframes(randsample(find(t.binLabels==x),sum(t.idxall)/t.binnum)), ...
                    unique(t.binLabels), 'UniformOutput', false));
            end
                conmat.mats.event_onset_frames(1,t.idxall) = t.poss_frames_mat(randperm(numel(t.poss_frames_mat)));
            
        end
    end
end

% check some troubling stuff
%max(conmat.mats.event_onset_frames(1,conmat.mats.eventnum==2))/p.scr_refrate;


% distribute second events randomly in possible interval
t.idx = find(~isnan(conmat.mats.eventRDK(2,:)));
for i_ev = 1:numel(t.idx)
    t.idxframes = find(t.poss_frames_2(conmat.mats.eventRDK(2,t.idx(i_ev)),:));
    % what is the time point of the first event
    t.idx2 = find(t.idxframes > conmat.mats.event_onset_frames(1,t.idx(i_ev))+p.scr_refrate*(p.stim.event.min_dist-0.01));
    conmat.mats.event_onset_frames(2,t.idx(i_ev)) = t.idxframes(t.idx2(randsample(numel(t.idx2),1)));
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

% include pre-cue events
% number of pre-cue events
t.mat = [repmat(p.stim.precue_event.num,1,floor(conmat.totaltrials/numel(p.stim.precue_event.num))), ...
    randsample(p.stim.precue_event.num,mod(conmat.totaltrials, numel(p.stim.precue_event.num)))];
conmat.mats.precue_eventnum = randsample(t.mat,numel(t.mat));

% indicate event: 1 - horizontal shorter, 2 - horizontal longer, 3 - vertical shorter, 4 - vertical longer
conmat.mats.precue_eventid = nan(size(conmat.mats.precue_eventnum));
t.mat = [repmat([1 2 3 4],1,floor(sum(conmat.mats.precue_eventnum)/4)),randsample([1 2 3 4],mod(sum(conmat.mats.precue_eventnum), 4))];
conmat.mats.precue_eventid(conmat.mats.precue_eventnum==1)= t.mat(randsample(numel(t.mat ),numel(t.mat )));

% target or distractor? 1 = target, 2= distractor
conmat.mats.precue_eventtype = nan(size(conmat.mats.precue_eventnum));
conmat.mats.precue_eventtype(ismember(conmat.mats.precue_eventid,p.stim.precue_event.targets))=1;
conmat.mats.precue_eventtype(ismember(conmat.mats.precue_eventid,setdiff([1 2 3 4], p.stim.precue_event.targets)))=2;

% pre-cue event onset times
t.idx = find(conmat.mats.precue_eventnum==1);
conmat.mats.precue_event_onset_times = nan(size(conmat.mats.precue_eventnum));
for i_tr = 1:numel(t.idx)
    % index possible positions
    t.idx1 = p.stim.precue_event.min_onset:0.1:conmat.mats.pre_cue_times(t.idx(i_tr))-(p.stim.precue_event.length+p.stim.precue_event.min_offset);
    conmat.mats.precue_event_onset_times(t.idx(i_tr)) = randsample(t.idx1,1);
end
conmat.mats.precue_event_onset_frames = round(conmat.mats.precue_event_onset_times*p.scr_refrate);

%% randomize all information across experiment
t.tidx = randperm(conmat.totaltrials);

t.fields = fieldnames(conmat.mats);
for i_fields = 1:numel(t.fields)
    conmat.mats.(t.fields{i_fields}) = conmat.mats.(t.fields{i_fields})(:,t.tidx);
end

conmat.mats.block = repmat(1:conmat.totalblocks,conmat.trialsperblock,1);
conmat.mats.block = conmat.mats.block(:)';

t.resp_hand =  {'linke','rechte'};
if flag_training == 1
    t.idx = randsample([1 2],1);
else
    t.idx = repmat([1 2],1,conmat.totalblocks/2);
end
conmat.mats.responsehand = repmat(t.resp_hand(t.idx(randperm(numel(t.idx)))),conmat.trialsperblock,1);
conmat.mats.responsehand = conmat.mats.responsehand(:)';

%% write all information into trial structure
% create frame mat, onset time for events
t.fields = fieldnames(conmat.mats);
conmat.trials = struct();

for i_field = 1:numel(t.fields)
    values = num2cell(conmat.mats.(t.fields{i_field}),[1]);
    [conmat.trials(1:conmat.totaltrials).(t.fields{i_field})] = values{:};              % Assign field dynamically
end
% add trial num
t.trials = num2cell(1:conmat.totaltrials);
[conmat.trials(1:conmat.totaltrials).trialnum] = t.trials{:};

% post-cue times
[conmat.trials(1:conmat.totaltrials).post_cue_times] = deal(p.stim.time_postcue);

% post-cue frames
[conmat.trials(1:conmat.totaltrials).post_cue_frames] = deal(p.stim.time_postcue*p.scr_refrate);

    

end

