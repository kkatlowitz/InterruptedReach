addpath(genpath('/Users/Dennis/Documents/MATLAB/Mogilner_Lab'));
addpath(genpath('/Users/Dennis/Documents/MATLAB/chronux_2_12'));
addpath(genpath('/Users/Dennis/Documents/MATLAB/wave_clus-testing'));

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

%% spike detection
par=set_parameters;
for s=1:1
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

%% this part finds the turn in the trajectories of turn and swerve trials in a supervised fashion

% you do not need to do this if you are loading previously selected times
% (see below)

% y: agree with automated choice (or previously selected choice if
%    scrolling back)
% n: disagree, select turn point with mouse and left-click
% q: there is no turn, void this trials turn time
% g: go back

lpfilt=designfilt('lowpassfir','FilterOrder',10,...  % low pass filter for trajectory to find turn
    'PassbandFrequency',1,'StopbandFrequency',25, ...
    'SampleRate',100);
lpfilt=lpfilt.Coefficients;

for s=[1,2,4,7]
    for blocknum=1:length(sub(s).block)
        b=sub(s).block(blocknum).b;
        turn_time=zeros(length(b.type),1);
        iter=1;
        turn_ind=zeros(1,length(sub(s).block(blocknum).b.t));
        while iter<=length(sub(s).block(blocknum).b.t)
            if b.type(iter)~=1 && ~isnan(b.event_time(iter,3)) %dont find the turn if the trial is a go trial or if the trial was a fixation break
                if turn_ind(iter)==0 || turn_ind(iter)==-1
                    [turn_time(iter),turn_ind(iter)]=knee(real(b.leap(iter).screen),... % find the turn
                        b.leap(iter).t,b.t(iter,1),lpfilt,[b.event_time(iter,3),inf],[200,300]);
                end
                plot(b.leap(iter).screen(:,1),b.leap(iter).screen(:,2));
                hold on
                scatter(b.leap(iter).screen(turn_ind(iter),1),b.leap(iter).screen(turn_ind(iter),2),...
                    'Marker','*','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',10)
                prompt=true;
                while prompt==true
                    %str=input('Is this the turn? (y/n/q)','s');
                    k=waitforbuttonpress;
                    if k==1
                        str=get(gcf,'CurrentCharacter');
                        if strcmp(str,'y')
                            prompt=false;
                            iter=iter+1;
                        elseif strcmp(str,'n')
                            prompt=false;
                            [x,y]=ginput(1);
                            xdiff=bsxfun(@minus,x,b.leap(iter).screen(:,1));
                            ydiff=bsxfun(@minus,y,b.leap(iter).screen(:,2));
                            [~,turn_ind(iter)]=min(sqrt(xdiff.^2+ydiff.^2));
                            turn_time(iter)=b.leap(iter).t(turn_ind(iter));
                            iter=iter+1;
                        elseif strcmp(str,'q')
                            prompt=false;
                            turn_time(iter)=-1;
                            turn_ind(iter)=-1;
                            iter=iter+1;
                        elseif strcmp(str,'g')
                            prompt=false;
                            prior_trials=1:iter;
                            iter=find(turn_ind(prior_trials)~=-1 & turn_ind(prior_trials)~=0,2,'last');
                            iter=iter(1);
                        else
                            continue
                        end
                    else
                        continue
                    end
                end
            else
                turn_time(iter)=NaN;
                iter=iter+1;
            end
            cla
        end
        close(gcf)
        turn_time_align=b.tev1_trials(:,3)+turn_time-b.event_time(:,3); %convert turn time to AO time
        turn_time_align(turn_ind==-1)=-1;
        b.tev1_trials=[b.tev1_trials,turn_time_align];
        sub(s).block(blocknum).b=b;
    end
end
%% extract the turn times that you so arduously selected along with the other 
%  relevant trial times into their own structure that can be saved and is
%  not massive (do not use if your are loading saved times)

for s=[1,2,4,7]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            sub_times(s).block(blocknum).b(lead).tev1_trials=sub(s).block(blocknum).b(lead).tev1_trials;
        end
    end
end

%% if you are loading previously saved times, run this block
load('sub_times')
for s=[1,2,4,7]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            sub(s).block(blocknum).b(lead).tev1_trials=sub_times(s).block(blocknum).b(lead).tev1_trials;
        end
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
for s=[1,2,4,7]
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
