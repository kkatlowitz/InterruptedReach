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

sub(8).pathbase=pathbase;
sub(8).pathbeh{1}='/Sub08/Block02/Sub08_WH_02052018_OR_HandR_STNL_Block02.mat';
sub(8).pathbeh{2}='/Sub08/Block03/Sub08_WH_02052018_OR_HandR_STNL_Block03.mat';
sub(8).pathneur{1}='/Sub08/Block02/LT1D4.315F0002.mat';
sub(8).pathneur{2}='/Sub08/Block03/LT1D2.894F0002.mat';

sub(9).pathbase=pathbase;
sub(9).pathbeh{1}='/Sub09/Block01/Sub09_PH_102418_OR_HandR_STNL_Block01.mat';
sub(9).pathneur{1}='/Sub09/Block01/LT1D15.000F000345.mat';
sub(9).pathbeh{2}='/Sub09/Block02/Sub09_PH_102418_OR_HandR_STNL_Block02.mat';
sub(9).pathneur{2}='/Sub09/Block02/LT1D2.062F0001.mat';

sub(10).pathbase=pathbase;
sub(10).pathbeh{1}='/Sub10/Block01/Sub10_SH_070119_OR_HandL_STNR_Block01.mat';
sub(10).pathneur{1}='/Sub10/Block01/RT1D4.671.mat';


sub(11).pathbase=pathbase;
sub(11).pathbeh{1}='/Sub11/Block01/Sub11_MO_070819_OR_HandR_STNL_Block1.mat';
sub(11).pathneur{1}='/Sub11/Block01/LT1D2.977.mat';
sub(11).pathbeh{2}='/Sub11/Block03/Sub11_MO_070819_OR_HandR_STNL_Block3.mat';
sub(11).pathneur{2}='/Sub11/Block03/LT1D1.267F0002.mat';
sub(11).pathbeh{3}='/Sub11/Block05/Sub11_MO_070819_OR_HandR_STNL_Block5.mat';
sub(11).pathneur{3}='/Sub11/Block05/LT1D0.251.mat';
sub(11).pathbeh{4}='/Sub11/Block06/Sub11_MO_080719_OR_HandL_STNR_Block6.mat';
sub(11).pathneur{4}='/Sub11/Block06/RT1D5.102.mat';
sub(11).pathbeh{5}='/Sub11/Block07/Sub11_MO_080719_OR_HandL_STNR_Block7.mat';
sub(11).pathneur{5}='/Sub11/Block07/RT1D1.506.mat';

%load deltafilt


tic
for s=[1,2,4,7,8,9,11]
    s
    for blocknum=1:length(sub(s).pathbeh)
        
        if s<7
            [ao,b]=extractData(pathbase,sub(s).pathneur{blocknum},sub(s).pathbeh{blocknum},[1,2],[],1);
        elseif s>=7 && s<10
            [ao,b]=extractData(pathbase,sub(s).pathneur{blocknum},sub(s).pathbeh{blocknum},[1],[],0);
        else
            [ao,b]=extractData(pathbase,sub(s).pathneur{blocknum},sub(s).pathbeh{blocknum},[1,2],[],0);
        end
        ao=filter_rawdata_v2(ao,pathbase,s,blocknum,SOS,G);
        ao.dat=[]; %get rid of the unfiltered data to save memory
        block.b=b;
        block.ao=ao;
        
        if blocknum==1
            %%sub(s).block(blocknum)=block;
        else
            %%sub(s).block=[sub(s).block,block];
        end
    end
    toc
end

%%sub(2).block(2).b.begin=sub(2).block(1).b.begin;
%sub(2).block(2).ao.lfp=[NaN*ones(length(sub(2).block(1).ao.lfp),2);sub(2).block(2).ao.lfp];

%this fixes an issue with the data file for Subject 2 Block 2
mfile1=matfile([pathbase,'/LFP_delta/Sub2_Block1_Lead1.mat']);
mfile2=matfile([pathbase,'/LFP_delta/Sub2_Block2_Lead1.mat']);
mfile2.Properties.Writable=true;
mfile2.lfp=[NaN*ones(size(mfile1,'lfp',1),1);mfile2.lfp];

mfile2=matfile([pathbase,'/LFP_delta/Sub2_Block2_Lead2.mat']);
mfile2.Properties.Writable=true;
mfile2.lfp=[NaN*ones(size(mfile1,'lfp',1),1);mfile2.lfp];

%% this part  finds the turn in the trajectories of turn and swerve trials in a supervised fashion

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
for s=11
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
                            turn_time(iter)=b.leap(iter).t(turn_ind(iter))-b.t(iter,1);
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

for s=[1,2,4,7,8,9,10,11]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            sub_times(s).block(blocknum).b(lead).tev1_trials=sub(s).block(blocknum).b(lead).tev1_trials;
        end
    end
end

%% select motion start times
lpfilt=designfilt('lowpassfir','FilterOrder',30,...  % low pass filter for trajectory to find turn
'PassbandFrequency',1,'StopbandFrequency',5, ...
'SampleRate',100);
%lpfilt=designfilt('lowpassiir','designmethod','butter','FilterOrder',30,...  % low pass filter for trajectory to find turn
%'HalfPowerFequency',5,'SampleRate',100);
lpfilt=lpfilt.Coefficients;

for s=11%[1,2,4,7,8,9,10,11]
    for blocknum=1:length(sub(s).block)
        b=sub(s).block(blocknum).b;
        mov_time=zeros(length(b.type),1);
        iter=1;
        mov_ind=zeros(1,length(sub(s).block(blocknum).b.t));
        while iter<=length(sub(s).block(blocknum).b.t)
            if ~isnan(b.event_time(iter,3))
                coords=b.leap(iter).screen;
                t=real(b.leap(iter).t);
                start_t=b.t(iter);
                if mov_ind(iter)==0 || mov_ind(iter)==-1
                    [mov_time(iter),mov_ind(iter)]=moving(coords,t,start_t,lpfilt,[b.event_time(iter,3),Inf]);
                end
                
                plot(b.leap(iter).t-b.t(iter,1),coords(:,1));
                hold on
                if ~isnan(mov_ind(iter))
                scatter(b.leap(iter).t(mov_ind(iter))-b.t(iter,1),b.leap(iter).screen(mov_ind(iter),1),...
                    'Marker','*','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',10)
                end
                title(num2str(iter))
                
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
                            xdiff=bsxfun(@minus,x,b.leap(iter).t-b.t(iter,1));
                            ydiff=bsxfun(@minus,y,b.leap(iter).screen(:,1));
                            xdiff=range(ydiff)/range(xdiff)*xdiff;
                            [~,mov_ind(iter)]=min(sqrt(xdiff.^2+ydiff.^2));
                            mov_time(iter)=b.leap(iter).t(mov_ind(iter))-b.t(iter,1);
                            iter=iter+1;
                        elseif strcmp(str,'q')
                            prompt=false;
                            mov_time(iter)=-1;
                            mov_ind(iter)=-1;
                            iter=iter+1;
                        elseif strcmp(str,'g')
                            prompt=false;
                            prior_trials=1:iter;
                            iter=find(mov_ind(prior_trials)~=-1 & mov_ind(prior_trials)~=0,2,'last');
                            iter=iter(1);
                        else
                            continue
                        end
                    else
                        continue
                    end
                end
            else
                mov_time(iter)=NaN;
                iter=iter+1;
            end
            cla
        end
        close(gcf)
        mov_time_align=b.tev1_trials(:,3)+mov_time-b.event_time(:,3); %convert move time to AO time
        mov_time_align(mov_ind==-1)=-1;
        b.tev1_trials=[b.tev1_trials,mov_time_align];
        sub(s).block(blocknum).b=b;
    end
end
%% if you are loading previously saved times, run this block
load('sub_times_v6')
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            sub(s).block(blocknum).b(lead).tev1_trials=sub_times(s).block(blocknum).b(lead).tev1_trials;
        end
    end
end
%% select artifacts
for s=[2]%[1,2,4,7,8,9]
    s
   for blocknum=2:length(sub(s).block)
       blocknum
               for lead=2%sub(s).block(blocknum).ao.leads
                   lead
                   mfile=matfile([pathbase,'/lfp/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(lead)]);
                   lfp=mfile.lfp;
                   artlog=artifact_sel_v2(lfp,sub(s).block(blocknum).b,2000);
                   sub(s).block(blocknum).ao.artlog(lead,:)=artlog;
                   %sub(s).block(blocknum).ao=artifact_sel(sub(s).block(blocknum).ao,sub(s).block(blocknum).b);
               end
   end
end
%% save artifacts
for s=2%[1,2,4,7,8,9,11]
    for blocknum=2%1:length(sub(s).block)
        arts(s).block(blocknum).artlog=sub(s).block(blocknum).ao.artlog;
    end
end
%% load previously saved artifacts
load arts_v6
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        sub(s).block(blocknum).ao.artlog=arts(s).block(blocknum).artlog;
    end
end

%% pseudo events
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        tev1_trials=sub(s).block(blocknum).b.tev1_trials;
        event_happened=sub(s).block(blocknum).b.event_happened;
        event_happened=[event_happened,tev1_trials(:,9)~=-1 & ~isnan(tev1_trials(:,9)),tev1_trials(:,10)~=-1 & ~isnan(tev1_trials(:,10))];
        event_happened(sub(s).block(blocknum).b.type==4,4)=0;
        event_happened(isnan(tev1_trials))=0;
        
        tev1_pseudo=tev1_trials;
        event_happened_pseudo=event_happened;
        pseudo_events=false(size(event_happened));
        
        trueevents=logical(event_happened(:,4));  %real turn signals
        ntrue=nnz(trueevents);
        duration=tev1_trials(:,8)-tev1_trials(:,10);
        truedist=(tev1_trials(trueevents,4)-tev1_trials(trueevents,10))./duration(trueevents); %distribution of real turn signals relative to movement-response duration
        
        falseevents=logical(event_happened(:,10)) & ~event_happened(:,4); %events where turn signal did not happen but movement did (i.e go swerve trials)
        nfalse=nnz(falseevents);
        tev1_pseudo(falseevents,4)=tev1_pseudo(falseevents,10)+truedist(unidrnd(ntrue,nfalse,1)).*duration(falseevents); %populate false event times from distribution of true event times
        event_happened_pseudo(falseevents,4)=1;
        pseudo_events(falseevents,4)=true;
        
        trueevents=logical(event_happened(:,9)); %real turns
        ntrue=nnz(trueevents);
        duration=tev1_trials(:,8)-tev1_pseudo(:,4);
        truedist=(tev1_pseudo(trueevents,9)-tev1_pseudo(trueevents,4))./duration(trueevents);
        
        falseevents=logical(event_happened_pseudo(:,4)) & ~event_happened(:,9); %events where there was a real or pseudo turn signal but no turn (i.e. all go trials and error turn/swerve trials)
        nfalse=nnz(falseevents);
        tev1_pseudo(falseevents,9)=tev1_pseudo(falseevents,4)+truedist(unidrnd(ntrue,nfalse,1)).*duration(falseevents); %populate false event times from distribution of true event times
        event_happened_pseudo(falseevents,9)=1;
        pseudo_events(falseevents,9)=true;
        
        sub(s).block(blocknum).bpseudo=sub(s).block(blocknum).b;
        sub(s).block(blocknum).bpseudo.tev1_trials=tev1_pseudo;
        sub(s).block(blocknum).bpseudo.event_happened=event_happened_pseudo;
        sub(s).block(blocknum).bpseudo.pseudo_events=pseudo_events;
    end
end
%% align the spikes and LFP to behavioral events
sub(1).pathbase='/Users/Dennis/Documents/MATLAB/Mogilner_Lab/interrupted_reach';

sub(1).pathbeh{1}='/Sub01/Block02/Sub01_DD_04192017_OR_HandL_STNR_Block02.mat';
sub(1).pathbeh{2}='/Sub01/Block04/Sub01_DD_04192017_OR_HandL_STNR_Block04.mat';

sub(1).pathneur{1}='/Sub01/Block02/RT1D4.152F0001.mat';
sub(1).pathneur{2}='/Sub01/Block04/RT1D2.682F0002.mat';

%sub(1).pathspikes{1,1}='/Sub01/Block02/times_Sub01_DD_04192017_OR_HandL_STNR_Block02_spikes1.mat';
sub(1).pathspikes{1,1}='/Sub01/Block02/plexon_Sub01_DD_04192017_OR_HandL_STNR_Block02_spikes1.txt';
%sub(1).pathspikes{1,2}='/Sub01/Block02/times_Sub01_DD_04192017_OR_HandL_STNR_Block02_spikes2.mat';
sub(1).pathspikes{1,2}='/Sub01/Block02/plexon_Sub01_DD_04192017_OR_HandL_STNR_Block02_spikes2.txt';
%sub(1).pathspikes{2,1}='/Sub01/Block04/times_Sub01_DD_04192017_OR_HandL_STNR_Block04_spikes1.mat';
sub(1).pathspikes{2,1}='/Sub01/Block04/plexon_Sub01_DD_04192017_OR_HandL_STNR_Block04_spikes1.txt';
%sub(1).pathspikes{2,2}='/Sub01/Block04/times_Sub01_DD_04192017_OR_HandL_STNR_Block04_spikes2.mat';
sub(1).pathspikes{2,2}='/Sub01/Block04/plexon_Sub01_DD_04192017_OR_HandL_STNR_Block04_spikes2.txt';

sub(2).pathbase=pathbase;
sub(2).pathbeh{1}='/Sub02/Block01/Sub02_AK_04262017_OR_HandR_STNL_Block01.mat';
sub(2).pathbeh{2}='/Sub02/Block02/Sub02_AK_04262017_OR_HandR_STNL_Block02.mat';

sub(2).pathneur{1}='/Sub02/Block01/LT1D5.074F0001_block1.mat';
sub(2).pathneur{2}='/Sub02/Block02/LT1D5.074F0001_block2.mat';

%sub(2).pathspikes{1,1}='/Sub02/Block01/times_Sub02_AK_04262017_OR_HandR_STNL_Block01_spikes1.mat';
sub(2).pathspikes{1,1}='/Sub02/Block01/plexon_Sub02_AK_04262017_OR_HandR_STNL_Block01_spikes1.txt';
%sub(2).pathspikes{2,1}='/Sub02/Block02/times_Sub02_AK_04262017_OR_HandR_STNL_Block02_spikes1.mat';
sub(2).pathspikes{2,1}='/Sub02/Block02/plexon_Sub02_AK_04262017_OR_HandR_STNL_Block02_spikes1.txt';
%sub(2).pathspikes{1,2}='/Sub02/Block01/times_Sub02_AK_04262017_OR_HandR_STNL_Block01_spikes2.mat';
sub(2).pathspikes{1,2}='/Sub02/Block01/plexon_Sub02_AK_04262017_OR_HandR_STNL_Block01_spikes2.txt';
%sub(2).pathspikes{2,2}='/Sub02/Block02/times_Sub02_AK_04262017_OR_HandR_STNL_Block02_spikes2.mat';
sub(2).pathspikes{2,2}='/Sub02/Block02/plexon_Sub02_AK_04262017_OR_HandR_STNL_Block02_spikes2.txt';

sub(4).pathbase=pathbase;
sub(4).pathbeh{1}='/Sub04/Block01/Sub04_EW_07242017_OR_HandL_STNR_Block01.mat';
sub(4).pathneur{1}='/Sub04/Block01/ao_data.mat';
%sub(4).pathspikes{1}='/Sub04/Block01/times_Sub04_EW_07242017_OR_HandL_STNR_Block01_spikes1.mat';
sub(4).pathspikes{1}='/Sub04/Block01/plexon_Sub04_EW_07242017_OR_HandL_STNR_Block01_spikes1.txt';

sub(7).pathbase=pathbase;
sub(7).pathbeh{1}='/Sub07/Block02/Sub07_JB_10162017_OR_HandR_STNL_Block02.mat';
sub(7).pathneur{1}='/Sub07/Block02/LT1D4.826F0001.mat';
%sub(7).pathspikes{1}='/Sub07/Block02/times_Sub07_JB_10162017_OR_HandR_STNL_Block02_spikes1.mat';
sub(7).pathspikes{1}='/Sub07/Block02/plexon_Sub07_JB_10162017_OR_HandR_STNL_Block02_spikes1.txt';

sub(8).pathbeh{1}='/Sub08/Block02/Sub08_WH_02052018_OR_HandR_STNL_Block02.mat';
sub(8).pathbeh{2}='/Sub08/Block03/Sub08_WH_02052018_OR_HandR_STNL_Block03.mat';

sub(8).pathneur{1}='/Sub08/Block02/LT1D4.315F0002.mat';
sub(8).pathneur{2}='/Sub08/Block03/LT1D2.894F0002.mat';

%sub(8).pathspikes{1,1}='/Sub08/Block02/times_Sub08_WH_02052018_OR_HandR_STNL_Block02_spikes1.mat';
sub(8).pathspikes{1,1}='/Sub08/Block02/plexon_Sub08_WH_02052018_OR_HandR_STNL_Block02_spikes1.txt';
%sub(8).pathspikes{2,1}='/Sub08/Block03/times_Sub08_WH_02052018_OR_HandR_STNL_Block03_spikes1.mat';
sub(8).pathspikes{2,1}='/Sub08/Block03/plexon_Sub08_WH_02052018_OR_HandR_STNL_Block03_spikes1.txt';


sub(9).pathbase=pathbase;
sub(9).pathbeh{1}='/Sub09/Block01/Sub09_PH_102418_OR_HandR_STNL_Block01.mat';
sub(9).pathbeh{2}='/Sub09/Block02/Sub09_PH_102418_OR_HandR_STNL_Block02.mat';

sub(9).pathneur{1}='/Sub09/Block01/LT1D15.000F000345.mat';
sub(9).pathneur{2}='/Sub09/Block02/LT1D2.062F0001.mat';

%sub(9).pathspikes{1,1}='/Sub09/Block02/times_Sub09_PH_102418_OR_HandR_STNL_Block02_spikes1.mat';
sub(9).pathspikes{1,1}='/Sub09/Block01/plexon_Sub09_PH_102418_OR_HandR_STNL_Block01_spikes1.txt';
sub(9).pathspikes{2,1}='/Sub09/Block02/plexon_Sub09_PH_102418_OR_HandR_STNL_Block02_spikes1.txt';

%sub(10).pathspikes{1,1}='/Sub10/Block01/plexon_Sub10_SH_07012019_OR_HandL_STNR_Block01_spikes1.txt';
%sub(10).pathspikes{1,2}='/Sub10/Block01/plexon_Sub10_SH_07012019_OR_HandL_STNR_Block01_spikes2.txt';

sub(10).pathbase=pathbase;
sub(10).pathbeh{1}='/Sub10/Block01/Sub10_SH_070119_OR_HandL_STNR_Block01.mat';
sub(10).pathneur{1}='/Sub10/Block01/RT1D4.671.mat';
sub(10).pathspikes{1,1}='/Sub10/Block01/plexon_Sub10_SH_070119_OR_HandL_STNR_Block01_spikes1.txt';
sub(10).pathspikes{1,2}='/Sub10/Block01/plexon_Sub10_SH_070119_OR_HandL_STNR_Block01_spikes2.txt';

sub(11).pathbase=pathbase;
sub(11).pathbeh{1}='/Sub11/Block01/Sub11_MO_070819_OR_HandR_STNL_Block1.mat';
sub(11).pathneur{1}='/Sub11/Block01/LT1D2.977.mat';
sub(11).pathspikes{1,1}='/Sub11/Block01/plexon_Sub11_MO_070819_OR_HandR_STNL_Block01_spikes1.txt';
sub(11).pathbeh{2}='/Sub11/Block03/Sub11_MO_070819_OR_HandR_STNL_Block3.mat';
sub(11).pathneur{2}='/Sub11/Block03/LT1D1.267F0002.mat';
sub(11).pathspikes{2,1}='/Sub11/Block03/plexon_Sub11_MO_070819_OR_HandR_STNL_Block03_spikes1.txt';
sub(11).pathspikes{2,2}='/Sub11/Block03/plexon_Sub11_MO_070819_OR_HandR_STNL_Block03_spikes2.txt';
sub(11).pathbeh{3}='/Sub11/Block05/Sub11_MO_070819_OR_HandR_STNL_Block5.mat';
sub(11).pathneur{3}='/Sub11/Block05/LT1D0.251.mat';
sub(11).pathspikes{3,1}='/Sub11/Block05/plexon_Sub11_MO_070819_OR_HandR_STNL_Block05_spikes1.txt';
sub(11).pathbeh{4}='/Sub11/Block06/Sub11_MO_080719_OR_HandL_STNR_Block6.mat';
sub(11).pathneur{4}='/Sub11/Block06/RT1D5.102.mat';
sub(11).pathspikes{4,1}='/Sub11/Block06/plexon_Sub11_MO_080719_OR_HandL_STNR_Block06_spikes1.txt';
sub(11).pathspikes{4,2}='/Sub11/Block06/plexon_Sub11_MO_080719_OR_HandL_STNR_Block06_spikes2.txt';
sub(11).pathbeh{5}='/Sub11/Block07/Sub11_MO_080719_OR_HandL_STNR_Block6.mat';
sub(11).pathneur{5}='/Sub11/Block07/RT1D1.506.mat';
sub(11).pathspikes{5,1}='/Sub11/Block07/plexon_Sub11_MO_080719_OR_HandL_STNR_Block07_spikes1.txt';
sub(11).pathspikes{5,2}='/Sub11/Block07/plexon_Sub11_MO_080719_OR_HandL_STNR_Block07_spikes2.txt';


clear b ao
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).pathbeh)
        
        for sp=1:size(sub(s).pathspikes,2)
            if isempty(sub(s).pathspikes{blocknum,sp})
                continue
            end
            cluster_class=importspikes([sub(s).pathbase,sub(s).pathspikes{blocknum,sp}]);
            b(sp)=align_spikes_v2(sub(s).block(blocknum).b(1),cluster_class);
            bpseudo(sp)=align_spikes_v2(sub(s).block(blocknum).bpseudo(1),cluster_class);
            
        end
        
        sub(s).block(blocknum).b=b;
        
        sub(s).block(blocknum).bpseudo=bpseudo;
        clear b bpseudo
    end
end
%% remove invalid cells and trials
for j=1:9
sub(4).block(1).b(1).su(j,1).data(sub(4).block(1).b(1).su(j,1).trials>49)=[];
sub(4).block(1).b(1).su(j,1).trials(sub(4).block(1).b(1).su(j,1).trials>49)=[];
end

for j=1:9
sub(11).block(1).b(1).su(j,3).data(sub(11).block(1).b(1).su(j,3).trials>27)=[];
sub(11).block(1).b(1).su(j,3).trials(sub(11).block(1).b(1).su(j,3).trials>27)=[];
end
sub(11).block(4).b(1).su(:,[2,3])=[];
for j=1:9
sub(11).block(5).b(2).su(j,3).data(sub(11).block(5).b(2).su(j,3).trials<51)=[];
sub(11).block(5).b(2).su(j,3).trials(sub(11).block(5).b(2).su(j,3).trials<51)=[];
end

for j=1:9
sub(4).block(1).bpseudo(1).su(j,1).data(sub(4).block(1).bpseudo(1).su(j,1).trials>49)=[];
sub(4).block(1).bpseudo(1).su(j,1).trials(sub(4).block(1).bpseudo(1).su(j,1).trials>49)=[];
end
for j=1:9
sub(11).block(1).bpseudo(1).su(j,3).data(sub(11).block(1).bpseudo(1).su(j,3).trials>27)=[];
sub(11).block(1).bpseudo(1).su(j,3).trials(sub(11).block(1).bpseudo(1).su(j,3).trials>27)=[];
end
sub(11).block(4).bpseudo(1).su(:,[2,3])=[];
for j=1:9
sub(11).block(5).bpseudo(2).su(j,3).data(sub(11).block(5).bpseudo(2).su(j,3).trials<51)=[];
sub(11).block(5).bpseudo(2).su(j,3).trials(sub(11).block(5).bpseudo(2).su(j,3).trials<51)=[];
end
%% import unit quality
sheetcount=1;
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).pathbeh)
        units=zeros(length(sub(s).block(blocknum).b),1);
        for lead=1:length(sub(s).block(blocknum).b)
            units(lead)=size(sub(s).block(blocknum).b(lead).su,2);
        end
        sheet=xlsread([pathbase,'/spikesorting.xlsx'],sheetcount,['A2:C',num2str(1+sum(units))]);
        for lead=1:length(sub(s).block(blocknum).b)
            unit=sheet(sheet(:,1)==lead,2)+1;
            quality=sheet(sheet(:,1)==lead,3);
            sub(s).block(blocknum).b(lead).quality=zeros(size(sub(s).block(blocknum).b(lead).su,2),1);
            sub(s).block(blocknum).b(lead).quality(unit)=quality;
        end
        sheetcount=sheetcount+1;
    end
end
%% CWT calc
tic
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            s
            blocknum
            lead
            
            
            data=load([pathbase,'/lfp/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(lead)]);
            lfp=data.lfp;
            wt=NaN*ones(55,length(lfp));
            
            artlog=sub(s).block(blocknum).ao.artlog(lead,:);
            if nnz(artlog)>1
                artinds=find(~artlog);
                artstart=artinds(find(diff(find(~artlog))>1))+1;
                artinds=find(artlog);
                artend=artinds(find(diff(find(artlog))>1));
                artend(end+1)=artinds(end);
                
                artstart=artstart-5000;
                artend=artend+5000;
                
                if artlog(1)==1
                    artstart=[1,artstart];
                end
                
                overlaps=find(artstart(2:end)<artend(1:end-1));
                artstart(overlaps+1)=[];
                artend(overlaps)=[];
                
                if artlog(1)==0
                    nonartstart=[1,artend+1];
                else
                    nonartstart=artend+1;
                    artstart(1)=[];
                end
                nonartend=[artstart-1,length(lfp)];
                
                for ii=1:length(nonartstart)
                    [wt(:,nonartstart(ii):nonartend(ii)),f]=cwt(lfp(nonartstart(ii):nonartend(ii)),2000,'NumOctaves',9,'VoicesPerOctave',6,'TimeBandwidth',20);
                end
            else
                [wt,f]=cwt(lfp,2000,'NumOctaves',9,'VoicesPerOctave',6,'TimeBandwidth',20);
            end

            wt=wt';
            save([pathbase,'/CWT/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(lead)],'wt','s','blocknum','lead','-v7.3');
            toc
        end
    end
end
%% Hilbert calc
tic
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            s
            blocknum
            lead
            
            
            data=load([pathbase,'/LFP_delta/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(lead)]);
            lfp=data.lfp;
            hil=NaN*ones(1,length(lfp));
            
            artlog=sub(s).block(blocknum).ao.artlog(lead,:);
            if nnz(artlog)>1
                artinds=find(~artlog);
                artstart=artinds(find(diff(find(~artlog))>1))+1;
                artinds=find(artlog);
                artend=artinds(find(diff(find(artlog))>1));
                artend(end+1)=artinds(end);
                
                artstart=artstart-5000;
                artend=artend+5000;
                
                if artlog(1)==1
                    artstart=[1,artstart];
                end
                
                overlaps=find(artstart(2:end)<artend(1:end-1));
                artstart(overlaps+1)=[];
                artend(overlaps)=[];
                
                if artlog(1)==0
                    nonartstart=[1,artend+1];
                else
                    nonartstart=artend+1;
                    artstart(1)=[];
                end
                nonartend=[artstart-1,length(lfp)];
                
                for ii=1:length(nonartstart)
                    [hil(1,nonartstart(ii):nonartend(ii))]=hilbert(lfp(nonartstart(ii):nonartend(ii)));
                end
            else
                hil(1,:)=hilbert(lfp);
            end

            hil=hil';
            save([pathbase,'/Hilbert/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(lead)],'hil','s','blocknum','lead','-v7.3');
            toc
        end
    end
end
