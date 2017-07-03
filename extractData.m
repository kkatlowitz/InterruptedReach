function [ao,b]=extractData(F,ao_file,behavioral_file,loadMicro,loadMacro)
%% Load alpha omega neural data
ao=struct();
ao.file=[F,ao_file];%store the data file for ref
if exist(ao.file,'file')
    variableInfo = who('-file', ao.file);
else
    error('File Not Found. Check yourself before you wreck yourself')
end
if loadMicro
    load([F,ao_file],'CRAW_01_KHz','CRAW_01_TimeBegin');%you need these for both
    if ismember('CRAW_01', variableInfo) 
        load([F,ao_file],'CRAW_01','CRAW_01_BitResolution','CRAW_01_Gain');
        lead1=double(CRAW_01)./(CRAW_01_BitResolution*CRAW_01_Gain);
    else
        lead1=[];
        disp('No Lead 1')
    end
    if ismember('CRAW_02', variableInfo) 
        load([F,ao_file],'CRAW_02','CRAW_02_BitResolution','CRAW_02_Gain');
        lead2=double(CRAW_02)./(CRAW_02_BitResolution*CRAW_02_Gain);
    else 
        disp('No Lead 2')
        lead2=[];
    end
        %convert raw data on micro electrodes
    ao.fs=CRAW_01_KHz*1000;%sampling rate
    ao.dat=[lead1;lead2];
else
    ao.fs=NaN;
    ao.dat=NaN;
    disp('skipped micro')
end
if loadMacro
    if ismember('CMacro_LFP_01', variableInfo)
        load([F,ao_file],'CMacro_LFP_01','CMacro_LFP_01_KHz','CMacro_LFP_01_BitResolution','CMacro_LFP_01_Gain');
        macro1=double(CMacro_LFP_01)./(CMacro_LFP_01_BitResolution*CMacro_LFP_01_Gain);%these numbers should be the same
        ao.lfp_fs=CMacro_LFP_01_KHz*1e3;
    else
        macro1=[];
        disp('No Macro 1')
    end
    if ismember('CMacro_LFP_02', variableInfo) 
        load([F,ao_file],'CMacro_LFP_02','CMacro_LFP_02_KHz','CMacro_LFP_02_BitResolution','CMacro_LFP_02_Gain');
        macro2=double(CMacro_LFP_02)./(CMacro_LFP_02_BitResolution*CMacro_LFP_02_Gain);
        ao.lfp_fs=CMacro_LFP_02_KHz*1e3;%should overwrite with same thing if we have both
    else
        macro2=[];
        disp('No Macro 2')
    end
    ao.macro=[macro1;macro2];
else
    ao.macro=NaN;
    ao.lfp_fs=NaN;
    disp('skipped macro')
end
%% Do Behavior


load([F,ao_file],'CPORT__1','CPORT__1_KHz');%load AO sync data
ev1=CPORT__1(2,:); % events recorded by AO
b.tev1=CPORT__1(1,:)/(CPORT__1_KHz*1000); % event timestamps
b.begin=CRAW_01_TimeBegin;
b.file=[F,behavioral_file];%store it
b.trial_start_ind=find(ev1==21002); %trial start times
b.trial_end_ind=find(ev1==21003); %trial end times
trial_id_ind=find(ev1==24874); %trial ID message, word 1

%in case we started recording in the middle of a trial, get rid of an end before a start
b.trial_end_ind(b.trial_end_ind<b.trial_start_ind(1))=[];
b.trial_start_ind(b.trial_start_ind>b.trial_end_ind(end))=[];
%if there is a start with no stop after it, remove it
nS=length(b.trial_start_ind);
rm_start=zeros(nS,1);
for i=1:nS-1
    currStart=b.trial_start_ind(i);
    nextStart=b.trial_start_ind(i+1);
    if ~sum(b.trial_end_ind>currStart&b.trial_end_ind<nextStart)
        rm_start(i)=1;
    end
end
b.trial_start_ind(logical(rm_start))=[];

%if there is a stop with no start after it, remove it
nS=length(b.trial_end_ind);
rm_end=zeros(nS,1);
for i=1:nS-1
    currEnd=b.trial_end_ind(i);
    nextEnd=b.trial_end_ind(i+1);
    if ~sum(b.trial_start_ind>currEnd&b.trial_start_ind<nextEnd)
        rm_end(i)=1;
    end
end
b.trial_end_ind(logical(rm_end))=[];

%get the index of the trial ID message
a1=bsxfun(@gt,trial_id_ind,b.trial_start_ind');%left side
a2=bsxfun(@lt,trial_id_ind,b.trial_end_ind');%right side
[trial_id_ind_ind,~]=find(a1 & a2); %find which trial each ID is in

possibleInds=[];
for i=1:length(b.trial_start_ind)
    possibleInds=[possibleInds,b.trial_start_ind(i):b.trial_end_ind(i)];
end
trial_id_ind=intersect(trial_id_ind,possibleInds);

%get the data message associated with each trial ID message
trial_num=zeros(1,length(trial_id_ind_ind));
trial_num(trial_id_ind_ind)=ev1(trial_id_ind(trial_id_ind_ind)+4);
% there is a glitch where when word 4
% (I still don't know what this word represents but it increases by 1
% each trial) is 20502, it is missed, so the data is shifted
glitch_trial_ind=find(ev1(trial_id_ind+3)==20501)+1;
trial_num(trial_id_ind_ind(glitch_trial_ind))=ev1(trial_id_ind(glitch_trial_ind)+3);
b.trial_id_ind_ind=trial_id_ind_ind;
%now we need to fix the broken bit in the trial numbering. this will be ad
%hoc for now
bitokend=find(trial_num==10111); %last trial number that's fine
trial_num(bitokend+1:end)=trial_num(bitokend+1:end)+128; %fix everything after that
b.trial_num=trial_num-10000;

% find the events
b.event_names={'fp_on','targ_on','fp_off','targ_chg','in_chg_window',...
    'in_knot1_window','response','feedback'};%for
% the
event_codes=[21004,21006,21005,21300,21301,21304,21012,21013];
b.events_ind_all=nan(length(trial_num),length(event_codes)); %trials by events
for i=1:length(event_codes)
    events_ind=find(ev1==event_codes(i)); %find the index of each event
%     events_ind(events_ind<b.trial_start_ind(1))=[]; %get rid of those before the start of the first recorded trial
%     events_ind(events_ind>b.trial_start_ind(end))=[];
    events_ind=intersect(events_ind,possibleInds);
    left=bsxfun(@gt,events_ind,b.trial_start_ind');
    right=bsxfun(@lt,events_ind,b.trial_end_ind');
    [events_ind_ind,~]=find(left & right); %find which trial each event is in
    try
        b.events_ind_all(events_ind_ind,i)=events_ind; %store the events
    catch
        disp('you have more than one event for each trial')
    end      
end
% now convert all event indices to event times
b.tev1_trials=NaN(length(trial_num),length(event_codes));
event_log=~isnan(b.events_ind_all);
% b.tev1_trials(event_log)=b.tev1(b.events_ind_all(event_log));
%% load the rome data
if ~isempty(behavioral_file)
    load([F,behavioral_file]);
    for trial=1:length(trial_data)
        b.rome_id(trial)=trial_data{trial}.trial_id;
        b.type(trial)=trial_data{trial}.trial_type;
        b.t(trial,1:2)=[trial_data{trial}.start_t,trial_data{trial}.end_t];
        b.leap(trial).t=trial_data{trial}.touch_data(:,4);
        b.leap(trial).hand_detected=trial_data{trial}.touch_data(:,3);
        b.leap(trial).screen=trial_data{trial}.touch_data(:,1:2);
        b.leap(trial).xyz=trial_data{trial}.touch_data(:,5:7);
        b.event_happened(trial,:)=trial_data{trial}.event_happened;
        b.event_time(trial,:)=trial_data{trial}.event_time;
    end
    r2ao=ismember(b.rome_id,b.trial_num);
    b.type=b.type(r2ao);
    b.t=b.t(r2ao,:);
    b.event_happened=b.event_happened(r2ao,:);
    b.event_time=b.event_time(r2ao,:);
    b.leap=b.leap(r2ao);
end