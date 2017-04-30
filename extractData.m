function [ao,b]=extractData(F,ao_file,behavioral_file,spikes_file,micro_leads,macro_leads)
%% Load alpha omega data
%too many variables, dont need to load all
ao=struct();
load([F,ao_file],'CPORT__1','CPORT__1_KHz');
ev1=CPORT__1(2,:); % events recorded by AO
tev1=CPORT__1(1,:)/(CPORT__1_KHz*1000); % event timestamps
load([F,ao_file],'CRAW_01','CRAW_02','CRAW_01_KHz',...
    'CRAW_01_BitResolution','CRAW_02_BitResolution','CRAW_01_Gain','CRAW_02_Gain','CRAW_01_TimeBegin');
% make bandpass filter

if ~isempty(micro_leads)
    bpFilt = designfilt('bandpassfir','FilterOrder',500, ...
        'CutoffFrequency1',250,'CutoffFrequency2',5000, ...
        'SampleRate',CRAW_01_KHz*1000);
    
    lead1=double(CRAW_01)./(CRAW_01_BitResolution*CRAW_01_Gain);
    lead2=double(CRAW_02)./(CRAW_02_BitResolution*CRAW_02_Gain);
    ao.file=[F,ao_file];
    ao.fs=CRAW_01_KHz*1000;
    ao.leads=micro_leads;
    ao.dat=[lead1;lead2]';
    ao.dat=ao.dat(:,micro_leads);
    ao.bp=filtfilt(bpFilt,ao.dat);
    %ao.lfp=resample(ao.dat,1,44,'spline');%switch this to a median filter
    %ao.fs_low=ao.fs/44;
end
if ~isempty(macro_leads)
    ao.cannula=1;
end


%%
load([F,behavioral_file]);
b.trial_start_ind=find(ev1==21002); %trial start times
b.trial_end_ind=find(ev1==21003); %trial end times
%in case we started recording in the middle of a trial, get rid of and end before a start
b.trial_end_ind(b.trial_end_ind<b.trial_start_ind(1))=[];
b.trial_start_ind(b.trial_start_ind>b.trial_end_ind(end))=[];


%I have no idea what this part is
b.trial_id_ind=find(ev1==24874); %trial ID message, word 1
[trial_id_ind_ind,~]=find(bsxfun(@gt,b.trial_id_ind,b.trial_start_ind') & bsxfun(@lt,b.trial_id_ind,b.trial_end_ind')); %find which trial each ID is in

trial_num=zeros(1,length(trial_id_ind_ind));
trial_num(trial_id_ind_ind)=ev1(b.trial_id_ind(trial_id_ind_ind)+4); %get the data message associated with each trial ID message

glitch_trial_ind=find(ev1(b.trial_id_ind+3)==20501)+1; % there is a glitch where when word 4 (I still don't know what this word represents but it increases by 1 each trial) is 20502, it is missed, so the data is shifted
trial_num(trial_id_ind_ind(glitch_trial_ind))=ev1(b.trial_id_ind(glitch_trial_ind)+3);
b.trial_id_ind_ind=trial_id_ind_ind;
%now we need to fix the broken bit in the trial numbering. this will be ad
%hoc for now

bitokend=find(trial_num==10111); %last trial number that's fine
trial_num(bitokend+1:end)=trial_num(bitokend+1:end)+128; %fix everything after that

b.trial_num=trial_num-10000;
%% find the events
b.event_names={'fp_on','targ_on','fp_off','targ_chg','in_chg_window','in_knot1_window','response','feedback'};%for
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

if ~isempty(spikes_file)
    
    load([F,spikes_file]);
    
    stimes=cluster_class(:,2)/1000 + CRAW_01_TimeBegin; %spike times converted to seconds
    
    b.psth_start_num=-2;
    b.psth_end_num=2;
    b.binwidth=0.1;
    
    b.trial_start=tev1(b.trial_start_ind)';
    b.trial_end=tev1(b.trial_end_ind)';
    
    event_plot=[1:6,8];
    for j=1:7
        b.psth_start=b.tev1_trials(:,event_plot(j)) + b.psth_start_num; %psth relative to fp_off
        b.psth_end=b.tev1_trials(:,event_plot(j))+ b.psth_end_num;
        
        trial_start_align=b.trial_start-b.psth_start + b.psth_start_num;
        trial_end_align=b.trial_end-b.psth_start + b.psth_start_num;
        
        binl=b.psth_start_num:0.002:(b.psth_end_num-b.binwidth); %sliding window
        binu=(b.psth_start_num+b.binwidth):0.002:b.psth_end_num;
        
        %ntrials=length(unique(su(1).strials));
        
        for i=1:length(unique(cluster_class(:,1)))-1
            b.su(j,i).stimes=stimes(cluster_class(:,1)==i); %single unit spike times
            
            b.su(j,i).stimes(b.su(j,i).stimes<tev1(b.trial_start_ind(1)))=[];
            
            [b.su(j,i).stimes_ind,b.su(j,i).strials]=find(bsxfun(@gt,b.su(j,i).stimes,b.psth_start')...
                & bsxfun(@le,b.su(j,i).stimes,b.psth_end') & bsxfun(@gt,b.su(j,i).stimes,b.trial_start')...
                & bsxfun(@le,b.su(j,i).stimes,b.trial_end'));
            
            b.su(j,i).stimes_align=b.su(j,i).stimes(b.su(j,i).stimes_ind)-b.psth_start(b.su(j,i).strials) + b.psth_start_num;
            b.su(j,i).trials=unique(b.su(j,i).strials);
            ntrials=length(unique(b.su(j,i).strials));
            
            
            spikebins=bsxfun(@ge,b.su(j,i).stimes_align,binl) & bsxfun(@lt,b.su(j,i).stimes_align,binu);
            spikebins=mat2cell(spikebins,histcounts(b.su(j,i).strials,[b.su(j,i).trials;b.su(j,i).strials(end)+0.5]));
            spikebins=cell2mat(cellfun(@(x)sum(x,1),spikebins,'UniformOutput',0));
            spikebins(~(bsxfun(@gt,binl,trial_start_align(b.su(j,i).trials)) & bsxfun(@le,binu,trial_end_align(b.su(j,i).trials))))=NaN;
            b.su(j,i).fr=mean(spikebins,'omitnan')/(b.binwidth);
            
            b.su(j,i).fr_se=std(spikebins,'omitnan')/(sqrt(length(b.su(j,i).trials))*b.binwidth);
            
            
        end
    end
end
