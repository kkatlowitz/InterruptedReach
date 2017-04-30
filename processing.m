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
        
        [ao,b]=extractData(pathbase,sub(s).pathneur{blocknum},sub(s).pathbeh{blocknum},[],1,[]);
        
        block=b;
        block.fs=ao.fs;
        block.bp=ao.bp;
        
        if blocknum==1
            sub(s).block(blocknum)=block;
        else
            sub(s).block=[sub(s).block,block];
        end
    end
end

%% spike detection
par=set_parameters;
for s=1:length(sub)
    for blocknum=1:length(sub(s).block)
        [spikes,index]=spike_adapt(sub(s).block(blocknum).bp,par);
        index=index/44;
        spikestruct.spikes=spikes;
        spikestruct.index=index;
        spikestruct.sr=44000;
        save([sub(s).pathbase,sub(s).pathbeh{blocknum}(1:end-4),'_spikes1.mat'],'-struct','spikestruct');
    end
end



%% now run wave_clus, save, and clear all

pathbase='/Users/Dennis/Documents/MATLAB/Mogilner_Lab/interrupted_reach';

sub(1).pathbeh{1}='/Sub01/Block02/Sub01_DD_04192017_OR_HandL_STNR_Block02.mat';
sub(1).pathbeh{2}='/Sub01/Block04/Sub01_DD_04192017_OR_HandL_STNR_Block04.mat';

sub(1).pathneur{1}='/Sub01/Block02/RT1D4.152F0001.mat';
sub(1).pathneur{2}='/Sub01/Block04/RT1D2.682F0002.mat';

sub(1).pathspikes{1}='/Sub01/Block02/times_Sub01_DD_04192017_OR_HandL_STNR_Block02_spikes1.mat';
sub(1).pathspikes{2}='/Sub01/Block04/times_Sub01_DD_04192017_OR_HandL_STNR_Block04_spikes1.mat';


for blocknum=1:length(sub(1).pathbeh)
    
    [~,b]=extractData(pathbase,sub(1).pathneur{blocknum},sub(1).pathbeh{blocknum},sub(1).pathspikes{blocknum},[],[]);
    
    sub(1).block(blocknum)=b;
    
end
