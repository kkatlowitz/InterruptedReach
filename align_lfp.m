function ao=align_lfp(ao,b)


b.psth_start_num=-4;
b.psth_end_num=4;
b.binwidth=0.1;

b.trial_start=b.tev1(b.trial_start_ind)';
b.trial_end=b.tev1(b.trial_end_ind)';

event_plot=[1:6,8,9];
ao.trials=cell(length(event_plot),size(b.tev1_trials,1));
for j=1:8
    psth_start=b.tev1_trials(:,event_plot(j)) + b.psth_start_num; %psth relative to fp_off
    psth_end=b.tev1_trials(:,event_plot(j))+ b.psth_end_num;
    
    trial_start_align=b.trial_start-psth_start + b.psth_start_num;
    trial_end_align=b.trial_end-psth_start + b.psth_start_num;
    
    psth_start_ind=mat2cell(round((psth_start-b.begin+1)*ao.fs_low),...
        ones(1,length(b.tev1_trials)));
    psth_end_ind=mat2cell(round((psth_end-b.begin+1)*ao.fs_low),...
        ones(1,length(b.tev1_trials)));
    
    padding=find(cellfun(@(x) x<=0, psth_start_ind));
    for i=1:length(padding)
        psth_start_ind{padding(i)}=1;
    end
    
    ao.trials(j,cellfun(@(x) ~isnan(x),psth_start_ind))=cellfun(@(x,y) ao.bp(x:y)',...
        psth_start_ind(cellfun(@(x) ~isnan(x),psth_start_ind)),...
        psth_end_ind(cellfun(@(x) ~isnan(x),psth_end_ind)),'UniformOutput',0);
    
    for i=1:length(padding)
        ao.trials{j,padding(i)}=padarray(ao.trials{j,padding(i)},[0,(b.psth_end_num-b.psth_start_num)*ao.fs_low+1-length(ao.trials{j,padding(i)})],0,'pre');
    end
    
end
end