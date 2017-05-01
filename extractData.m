function [ao,b]=extractData(F,ao_file,behavioral_file,micro_leads,macro_leads)
%% Load alpha omega neural data
ao=struct();
ao.file=[F,ao_file];%store the data file for ref
load([F,ao_file],'CRAW_01','CRAW_02','CRAW_01_KHz',...
    'CRAW_01_BitResolution','CRAW_02_BitResolution','CRAW_01_Gain','CRAW_02_Gain','CRAW_01_TimeBegin');
if ~isempty(micro_leads)
    %convert raw data on micro electrodes
    lead1=double(CRAW_01)./(CRAW_01_BitResolution*CRAW_01_Gain);
    lead2=double(CRAW_02)./(CRAW_02_BitResolution*CRAW_02_Gain);
    ao.fs=CRAW_01_KHz*1000;%sampling rate
    ao.leads=micro_leads;
    ao.dat=[lead1;lead2]';
    ao.dat=ao.dat(:,micro_leads);
end
if ~isempty(macro_leads)
    ao.cannula=1;
    %need to look up which is the LFP sig
end
%% Do Behavior


% DENNIS: WHAT OF THIS DO WE ACTUALLY NEED?
load([F,ao_file],'CPORT__1','CPORT__1_KHz');%load AO sync data
ev1=CPORT__1(2,:); % events recorded by AO
tev1=CPORT__1(1,:)/(CPORT__1_KHz*1000); % event timestamps
b.file=[F,behavioral_file];%store it
b.trial_start_ind=find(ev1==21002); %trial start times
b.trial_end_ind=find(ev1==21003); %trial end times
%in case we started recording in the middle of a trial, get rid of and end before a start
b.trial_end_ind(b.trial_end_ind<b.trial_start_ind(1))=[];
b.trial_start_ind(b.trial_start_ind>b.trial_end_ind(end))=[];


%get the index of the ?
trial_id_ind=find(ev1==24874); %trial ID message, word 1
[trial_id_ind_ind,~]=find(bsxfun(@gt,trial_id_ind,b.trial_start_ind') & bsxfun(@lt,trial_id_ind,b.trial_end_ind')); %find which trial each ID is in

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
    events_ind(events_ind<b.trial_start_ind(1))=[]; %get rid of those before the start of the first recorded trial
    events_ind(events_ind>b.trial_start_ind(end))=[];
    [events_ind_ind,~]=find(bsxfun(@gt,events_ind,b.trial_start_ind') & bsxfun(@lt,events_ind,b.trial_end_ind')); %find which trial each event is in
    b.events_ind_all(events_ind_ind,i)=events_ind; %store the events
end
% now convert all event indices to event times
b.tev1_trials=NaN(length(trial_num),length(event_codes));
event_log=~isnan(b.events_ind_all);
b.tev1_trials(event_log)=tev1(b.events_ind_all(event_log));
%% load the rome data
load([F,behavioral_file]);
for trial=1:length(trial_data)
    b.id(trial)=trial_data{trial}.id;
    b.t(trial,1:2)=[trial_data{trial}.start_t,trial_data{trial}.end_t];
    b.leap(trial).t=trial_data{trial}.touch_data(:,4);
    b.leap(trial).hand_detected=trial_data{trial}.touch_data(:,3);
    b.leap(trial).screen=trial_data{trial}.touch_data(:,1:2);
    b.leap(trial).xyz=trial_data{trial}.touch_data(:,5:7);
    b.event_happened(trial,:)=trial_data{trial}.event_happened;
    b.event_time(trial,:)=trial_data{trial}.event_time;
end
