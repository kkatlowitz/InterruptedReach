load('/Users/Dennis/Documents/MATLAB/Mogilner_Lab/interrupted_reach/Sub01/Block04/RT1D2.682F0002.mat')
load('/Users/Dennis/Documents/MATLAB/Mogilner_Lab/interrupted_reach/Sub01/Block04/Sub01_DD_04192017_OR_HandL_STNR_Block04.mat')
%%
ev1=CPORT__1(2,:); % events recorded by AO
tev1=CPORT__1(1,:)/(CPORT__1_KHz*1000); % event timestamps

t=(0:length(CRAW_01)-1)/(CRAW_01_KHz*1000);
lead1=double(CRAW_01)./(CRAW_01_BitResolution*CRAW_01_Gain);
lead2=double(CRAW_02)./(CRAW_02_BitResolution*CRAW_02_Gain);
%%
%lead1=filtfilt(bp300_3000,1,lead1);
block2.data=lead1;
block2.sr=CRAW_01_KHz*1000;
save('/Users/Dennis/Documents/MATLAB/Mogilner_Lab/interrupted_reach/Sub01/Block04/Sub01_DD_04192017_OR_HandL_STNR_Block04_MUR1.mat','-struct','block2');
block2.data=lead2;
block2.sr=CRAW_02_KHz*1000;
save('/Users/Dennis/Documents/MATLAB/Mogilner_Lab/interrupted_reach/Sub01/Block04/Sub01_DD_04192017_OR_HandL_STNR_Block04_MUR2.mat','-struct','block2');
% now run wave_clus
%%
load('/Users/Dennis/Documents/MATLAB/Mogilner_Lab/interrupted_reach/Sub01/Block04/times_Sub01_DD_04192017_OR_HandL_STNR_Block04_MUR1.mat')
%% find the trials
trial_start_ind=find(ev1==21002); %trial start times
trial_end_ind=find(ev1==21003); %trial end times
trial_end_ind(trial_end_ind<trial_start_ind(1))=[];%in case we started recording in the middle of a trial, get rid of and end before a start
trial_start_ind(trial_start_ind>trial_end_ind(end))=[];

trial_id_ind=find(ev1==24874); %trial ID message, word 1
[trial_id_ind_ind,~]=find(bsxfun(@gt,trial_id_ind,trial_start_ind') & bsxfun(@lt,trial_id_ind,trial_end_ind')); %find which trial each ID is in

trial_num=zeros(1,length(trial_id_ind_ind));
trial_num(trial_id_ind_ind)=ev1(trial_id_ind(trial_id_ind_ind)+4); %get the data message associated with each trial ID message

glitch_trial_ind=find(ev1(trial_id_ind+3)==20501)+1; % there is a glitch where when word 4 (I still don't know what this word represents but it increases by 1 each trial) is 20502, it is missed, so the data is shifted
trial_num(trial_id_ind_ind(glitch_trial_ind))=ev1(trial_id_ind(glitch_trial_ind)+3);

%now we need to fix the broken bit in the trial numbering. this will be ad
%hoc for now

bitokend=find(trial_num==10111); %last trial number that's fine
trial_num(bitokend+1:end)=trial_num(bitokend+1:end)+128; %fix everything after that

trial_num=trial_num-10000;
%% find the events
event_names={'fp_on','targ_on','fp_off','targ_chg','in_chg_window','in_knot1_window','response','feedback'};

event_codes=[21004,21006,21005,21300,21301,21304,21012,21013];

events_ind_all=NaN.*ones(length(trial_num),length(event_codes)); %trials by events
for i=1:length(event_codes)
    events_ind=find(ev1==event_codes(i)); %find the index of each event
    events_ind(events_ind<trial_start_ind(1))=[]; %get rid of those before the start of the first recorded trial
    events_ind(events_ind>trial_start_ind(end))=[];
    [events_ind_ind,~]=find(bsxfun(@gt,events_ind,trial_start_ind') & bsxfun(@lt,events_ind,trial_end_ind')); %find which trial each event is in
    events_ind_all(events_ind_ind,i)=events_ind; %store the events
end

% now convert all event indices to event times
tev1_trials=NaN.*ones(length(trial_num),length(event_codes));
event_log=~isnan(events_ind_all);
tev1_trials(event_log)=tev1(events_ind_all(event_log));
%% lets make some rasters

colors=get(gca,'colororder');
stimes=cluster_class(:,2)/1000 + CRAW_01_TimeBegin; %spike times converted to seconds

psth_start_num=-1;
psth_end_num=1;
psth_start=tev1_trials(:,3)-1; %psth relative to fp_off
psth_end=tev1_trials(:,3)+1;

trial_start=tev1(trial_start_ind)';
trial_end=tev1(trial_end_ind)';

trial_start_align=trial_start-psth_start-1;
trial_end_align=trial_end-psth_start-1;

binwidth=0.05;
binl=psth_start_num:0.002:(psth_end_num-binwidth); %sliding window
binu=(psth_start_num+binwidth):0.002:psth_end_num;

%ntrials=length(unique(su(1).strials));

for i=1:3
    su(i).stimes=stimes(cluster_class(:,1)==i); %single unit spike times
    
    su(i).stimes(su(i).stimes<tev1(trial_start_ind(1)))=[];
    
    [su(i).stimes_ind,su(i).strials]=find(bsxfun(@gt,su(i).stimes,psth_start')...
        & bsxfun(@le,su(i).stimes,psth_end') & bsxfun(@gt,su(i).stimes,trial_start')...
        & bsxfun(@le,su(i).stimes,trial_end'));
    
    su(i).stimes_align=su(i).stimes(su(i).stimes_ind)-psth_start(su(i).strials)-1;
    su(i).trials=unique(su(i).strials);
    ntrials=length(unique(su(1).strials));
    subplot(2,3,i)
    plot(su(i).stimes_align*1000,su(i).strials,'.') %plot rasters
    title(['Single Unit ',num2str(i)])
    xlabel('Time (ms)')
    set(gca,'box','off')
    subplot(2,3,3+i)
    spikebins=bsxfun(@ge,su(i).stimes_align,binl) & bsxfun(@lt,su(i).stimes_align,binu);
    spikebins=mat2cell(spikebins,histcounts(su(i).strials,[su(i).trials;su(i).strials(end)+0.5]));
    spikebins=cell2mat(cellfun(@(x)sum(x,1),spikebins,'UniformOutput',0));
    spikebins(~(bsxfun(@gt,binl,trial_start_align(su(i).trials)) & bsxfun(@le,binu,trial_end_align(su(i).trials))))=NaN;
    su(i).fr=mean(spikebins,'omitnan')/(binwidth);
    %su(i).fr=sum(bsxfun(@ge,su(i).stimes_align,binl) & bsxfun(@lt,su(i).stimes_align,binu))/(ntrials*binwidth);
    
    su(i).fr_se=std(spikebins,'omitnan')/(sqrt(length(su(i).trials))*binwidth);
    
    patchx=[(binl+binu)/2*1000,fliplr((binl+binu)/2*1000)];
    patchy=[su(i).fr-su(i).fr_se,fliplr(su(i).fr+su(i).fr_se)];
    patch(patchx(~isnan(patchy)),patchy(~isnan(patchy)),mean([colors(1,:);1,1,1]),'EdgeColor','none')
    hold on
    plot((binl+binu)/2*1000,su(i).fr)
    line([0,0],[0,100],'color','r')
    ylim([0,100])
    xlabel('Time (ms)')
    ylabel('Firing Rate (Hz)')
    set(gca,'box','off')
end

%% reformatting the behavioral data

result=cell(length(trial_data),1);
type=result;
for i=1:length(trial_data)
    result(i)=trial_data{i}.result;
    type(i)=trial_data{i}.trial_type;
end

