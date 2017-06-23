pathbase='/Users/Dennis/Documents/MATLAB/Mogilner_Lab/interrupted_reach';
sub(1).pathbase=pathbase;
sub(1).pathbeh{1}='/Sub01/Block02/Sub01_DD_04192017_OR_HandL_STNR_Block02.mat';
sub(1).pathbeh{2}='/Sub01/Block04/Sub01_DD_04192017_OR_HandL_STNR_Block04.mat';

sub(1).pathneur{1}='/Sub01/Block02/RT1D4.152F0001.mat';
sub(1).pathneur{2}='/Sub01/Block04/RT1D2.682F0002.mat';

%sub(2).pathbeh{1}='/Sub02/Block01/Sub02_AK_04262017_OR_HandR_STNL_Block01.mat';
%sub(2).pathbeh{1}='/Sub02/Block02/Sub02_AK_04262017_OR_HandR_STNL_Block02.mat';

%sub(2).pathneur{1}='/Sub02/Block01/LT1D14.009F0001.mat';
%sub(2).pathneur{1}='/Sub02/Block02/LT1D5.074F0002.mat';

for s=1:length(sub)
    for blocknum=1:length(sub(s).pathbeh)
        
        [ao,b]=extractData(pathbase,sub(s).pathneur{blocknum},sub(s).pathbeh{blocknum},1,[]);
        ao=filter_rawdata(ao);
        ao.dat=[]; %get rid of the unfiltered data to save memory
        block.b=b;
        block.ao=ao;
        
        if blocknum==1
            sub(s).block(blocknum)=block;
        else
            sub(s).block=[sub(s).block,block];
        end
    end
end

%% spike detection: only run this the first time and save the results with wave_clus. The next section will then just pull from the saved spikes
par=set_parameters;
for s=1:length(sub)
    for blocknum=1:length(sub(s).block)
        [spikes,index]=spike_adapt(sub(s).block(blocknum).ao.bp,par);
        index=index/44;
        spikestruct.spikes=spikes;
        spikestruct.index=index;
        spikestruct.sr=44000;
        save([sub(s).pathbase,sub(s).pathbeh{blocknum}(1:end-4),'_spikes1.mat'],'-struct','spikestruct');
    end
end



% if you are spike detecting, run wave_clus, save, clear all, and run this script from the beginning omitting the spike detection step.

%% detect turning
lpfilt=designfilt('lowpassfir','FilterOrder',10,...  % low pass filter for trajectory to find turn
    'PassbandFrequency',1,'StopbandFrequency',25, ...
    'SampleRate',100);
lpfilt=lpfilt.Coefficients;

for s=1:length(sub)
    for blocknum=1:length(sub(s).block)
        b=sub(s).block(blocknum).b;
        turn_time=zeros(length(b.type),1);
        for i=1:length(sub(s).block(blocknum).b.t)
            if b.type(i)~=1 && ~isnan(b.event_time(i,3)) %dont find the turn if the trial is a go trial or if the trial was a fixation break
                turn_time(i)=knee(real(b.leap(i).screen),... % find the turn
                    b.leap(i).t,b.t(i,1),lpfilt,[b.event_time(i,3),inf],[200,300]); 
            else
                turn_time(i)=NaN;
            end
        end
        turn_time_align=b.tev1_trials(:,3)+turn_time-b.event_time(:,3); %covnert turn time to AO time
        b.tev1_trials=[b.tev1_trials,turn_time_align];
        sub(s).block(blocknum).b=b;
    end
end


%% align the spikes and LFP to behavioral events
sub(1).pathbase='/Users/Dennis/Documents/MATLAB/Mogilner_Lab/interrupted_reach';

sub(1).pathbeh{1}='/Sub01/Block02/Sub01_DD_04192017_OR_HandL_STNR_Block02.mat';
sub(1).pathbeh{2}='/Sub01/Block04/Sub01_DD_04192017_OR_HandL_STNR_Block04.mat';

sub(1).pathneur{1}='/Sub01/Block02/RT1D4.152F0001.mat';
sub(1).pathneur{2}='/Sub01/Block04/RT1D2.682F0002.mat';

sub(1).pathspikes{1}='/Sub01/Block02/times_Sub01_DD_04192017_OR_HandL_STNR_Block02_spikes1.mat';
sub(1).pathspikes{2}='/Sub01/Block04/times_Sub01_DD_04192017_OR_HandL_STNR_Block04_spikes1.mat';



for blocknum=1:length(sub(1).pathbeh)
    
    load([sub(1).pathbase,sub(1).pathspikes{blocknum}]);
    
    b=align_spikes(sub(1).block(blocknum).b,cluster_class);
    ao=align_lfp(sub(1).block(blocknum).ao,sub(1).block(1).b);
    
    sub(1).block(blocknum).b=b;
    sub(1).block(blocknum).ao=ao;
    
end
