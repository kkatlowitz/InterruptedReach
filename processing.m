pathbase='/Users/Dennis/Documents/MATLAB/Mogilner_Lab/interrupted_reach';
sub(1).pathbase=pathbase;
sub(1).pathbeh{1}='/Sub01/Block02/Sub01_DD_04192017_OR_HandL_STNR_Block02.mat';
sub(1).pathbeh{2}='/Sub01/Block04/Sub01_DD_04192017_OR_HandL_STNR_Block04.mat';

sub(1).pathneur{1}='/Sub01/Block02/RT1D4.152F0001.mat';
sub(1).pathneur{2}='/Sub01/Block04/RT1D2.682F0002.mat';

sub(2).pathbase=pathbase;
sub(2).pathbeh{1}='/Sub02/Block01/Sub02_AK_04262017_OR_HandR_STNL_Block01.mat';
sub(2).pathbeh{2}='/Sub02/Block02/Sub02_AK_04262017_OR_HandR_STNL_Block02.mat';

sub(2).pathneur{1}='/Sub02/Block01/LT1D5.074F0001_block1.mat';
sub(2).pathneur{2}='/Sub02/Block02/LT1D5.074F0001_block2.mat';

sub(4).pathbase=pathbase;
sub(4).pathbeh{1}='/Sub04/Block01/Sub04_EW_07242017_OR_HandL_STNR_Block01.mat';
sub(4).pathneur{1}='/Sub04/Block01/ao_data.mat';

sub(7).pathbase=pathbase;
sub(7).pathbeh{1}='/Sub07/Block02/Sub07_JB_10162017_OR_HandR_STNL_Block02.mat';
sub(7).pathneur{1}='/Sub07/Block02/LT1D4.826F0001.mat';


for s=1:length(sub)
    for blocknum=1:length(sub(s).pathbeh)
        
        if s~=7
            [ao,b]=extractData(pathbase,sub(s).pathneur{blocknum},sub(s).pathbeh{blocknum},[1,2],[]);
        else
            [ao,b]=extractData(pathbase,sub(s).pathneur{blocknum},sub(s).pathbeh{blocknum},[1],[]);
        end
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

%% spike detection: no need to do if you already have spikes saved
par=set_parameters;
for s=1:7
    for blocknum=1:length(sub(s).block)
        [spikes,index]=spike_adapt(sub(s).block(blocknum).ao.bp(:,1),par);
        index=index/44;
        spikestruct.spikes=spikes;
        spikestruct.index=index;
        spikestruct.sr=44000;
        save([sub(s).pathbase,sub(s).pathbeh{blocknum}(1:end-4),'_spikes1.mat'],'-struct','spikestruct');
    end
end



%% now run wave_clus, save, and clear all

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

sub(1).pathspikes{1,1}='/Sub01/Block02/times_Sub01_DD_04192017_OR_HandL_STNR_Block02_spikes1.mat';
sub(1).pathspikes{1,2}='/Sub01/Block02/times_Sub01_DD_04192017_OR_HandL_STNR_Block02_spikes2.mat';
sub(1).pathspikes{2,1}='/Sub01/Block04/times_Sub01_DD_04192017_OR_HandL_STNR_Block04_spikes1.mat';
sub(1).pathspikes{2,2}='/Sub01/Block04/times_Sub01_DD_04192017_OR_HandL_STNR_Block04_spikes2.mat';

sub(2).pathbase=pathbase;
sub(2).pathbeh{1}='/Sub02/Block01/Sub02_AK_04262017_OR_HandR_STNL_Block01.mat';
sub(2).pathbeh{2}='/Sub02/Block02/Sub02_AK_04262017_OR_HandR_STNL_Block02.mat';

sub(2).pathneur{1}='/Sub02/Block01/LT1D5.074F0001_block1.mat';
sub(2).pathneur{2}='/Sub02/Block02/LT1D5.074F0001_block2.mat';

sub(2).pathspikes{1,1}='/Sub02/Block01/times_Sub02_AK_04262017_OR_HandR_STNL_Block01_spikes1.mat';
sub(2).pathspikes{2,1}='/Sub02/Block02/times_Sub02_AK_04262017_OR_HandR_STNL_Block02_spikes1.mat';
sub(2).pathspikes{1,2}='/Sub02/Block01/times_Sub02_AK_04262017_OR_HandR_STNL_Block01_spikes2.mat';
sub(2).pathspikes{2,2}='/Sub02/Block02/times_Sub02_AK_04262017_OR_HandR_STNL_Block02_spikes2.mat';

sub(4).pathbase=pathbase;
sub(4).pathbeh{1}='/Sub04/Block01/Sub04_EW_07242017_OR_HandL_STNR_Block01.mat';
sub(4).pathneur{1}='/Sub04/Block01/ao_data.mat';
sub(4).pathspikes{1}='/Sub04/Block01/times_Sub04_EW_07242017_OR_HandL_STNR_Block01_spikes1.mat';

sub(7).pathbase=pathbase;
sub(7).pathbeh{1}='/Sub07/Block02/Sub07_JB_10162017_OR_HandR_STNL_Block02.mat';
sub(7).pathneur{1}='/Sub07/Block02/LT1D4.826F0001.mat';
sub(7).pathspikes{1}='/Sub07/Block02/times_Sub07_JB_10162017_OR_HandR_STNL_Block02_spikes1.mat';

clear b ao
for s=1:length(sub)
    for blocknum=1:length(sub(s).pathbeh)
        
        for sp=1:size(sub(s).pathspikes,2)
            load([sub(s).pathbase,sub(s).pathspikes{blocknum,sp}]);
            
            b(sp)=align_spikes(sub(s).block(blocknum).b(1),cluster_class);
            ao_init=sub(s).block(blocknum).ao;
            ao_init.lfp=ao_init.lfp(:,sp);
            ao(sp)=align_lfp(ao_init,sub(s).block(blocknum).b);
            
        end
        
        sub(s).block(blocknum).b=b;
        sub(s).block(blocknum).ao=ao;
        clear ao b
    end
end
