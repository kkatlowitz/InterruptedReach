function b=align_spikes(b,cluster_class)

stimes=cluster_class(:,2)/1000 + b.begin; %spike times converted to seconds

b.psth_start_num=-2;
b.psth_end_num=2;
b.binwidth=0.1;

b.trial_start=b.tev1(b.trial_start_ind)';
b.trial_end=b.tev1(b.trial_end_ind)';

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

        b.su(j,i).stimes(b.su(j,i).stimes<b.tev1(b.trial_start_ind(1)))=[];

        [b.su(j,i).stimes_ind,b.su(j,i).strials]=find(bsxfun(@gt,b.su(j,i).stimes,b.psth_start')...
            & bsxfun(@le,b.su(j,i).stimes,b.psth_end') & bsxfun(@gt,b.su(j,i).stimes,b.trial_start')...
            & bsxfun(@le,b.su(j,i).stimes,b.trial_end'));

        b.su(j,i).stimes_align=b.su(j,i).stimes(b.su(j,i).stimes_ind)-b.psth_start(b.su(j,i).strials) + b.psth_start_num;
        b.su(j,i).trials=unique(b.su(j,i).strials);
        
        b.su(j,i).data=cell2struct(mat2cell(b.su(j,i).stimes_align,histcounts(b.su(j,i).strials,[b.su(j,i).trials;b.su(j,i).strials(end)+0.5])),'times',2);
    end
end
