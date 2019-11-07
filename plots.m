%% calculate PSTHs separated by trial type (requires chronux toolbox)
event_ord=[1,2,3,10,4,9,8];
suevent_ord=[1,2,3,9,4,8,7];
types=[1,4,5];

preadj=[0.5,0.05,0.05,0.05,0.05,0.05,0.05];
postadj=[0,0,-0.1,0,-0.1,0,0];

preall=cell(11,length(event_ord));
postall=preall;
typeallall=preall;
fixbreakallall=preall;


presplit=cell(11,5,2,length(event_ord),3);
postsplit=presplit;
typeall=presplit;
fixbreakall=presplit;
tic
for s=[1,2,4,7,8,9,11]
    s
    for blocknum=1:length(sub(s).block)
        blocknum
        result_orig=sub(s).block(blocknum).b(1).result;
        type_orig=sub(s).block(blocknum).b(1).type;
        trialtimes_orig=sub(s).block(blocknum).b(1).t;
        tev1_trials_orig=sub(s).block(blocknum).b(1).tev1_trials;
        tev1_trials_orig(type_orig==4,4)=NaN;
        event_happened_orig=sub(s).block(blocknum).b(1).event_happened;
        event_happened_orig=[event_happened_orig,tev1_trials_orig(:,9)~=-1 & ~isnan(tev1_trials_orig(:,9)),tev1_trials_orig(:,10)~=-1 & ~isnan(tev1_trials_orig(:,10))];
        event_happened_orig(type_orig==4,4)=0;
        event_happened_orig=event_happened_orig(:,event_ord);
        tev1_trials_orig=tev1_trials_orig(:,event_ord);
        %result_orig(result_orig==2 & type_orig==4 & ~(tev1_trials_orig(:,6)==-1 | isnan(tev1_trials_orig(:,6)))')=1; %turn trials where a turn was identified but were marked incorrect
        %result_orig(result_orig==2 & type_orig==5 & ~isnan(tev1_trials_orig(:,5))' & ~(tev1_trials_orig(:,6)==-1 | isnan(tev1_trials_orig(:,6)))')=1; %swerve trials where a swerve signal was present and a turn was identified but were marked incorrect
        result_orig(result_orig==1 & (type_orig==4 | type_orig==5) & tev1_trials_orig(:,6)'==-1)=2; %turn or swerve trials marked correct but no turn identified
        result_orig(~isnan(tev1_trials_orig(:,6)) & tev1_trials_orig(:,6)~=-1 & tev1_trials_orig(:,6)<tev1_trials_orig(:,3))=4; %trials where a turn was identified before the go cue (i.e. incorrectly identified), ideally would go back and fix turn identification here
        tev1_trials_orig(result_orig'==1 & tev1_trials_orig(:,4)<tev1_trials_orig(:,3),4)=tev1_trials_orig(result_orig'==1 & tev1_trials_orig(:,4)<tev1_trials_orig(:,3),3)+0.01;
        %fixbreak_orig=tev1_trials_orig(:,4)<tev1_trials_orig(:,3) | isnan(tev1_trials_orig(:,2)) | isnan(tev1_trials_orig(:,3));
        %error_orig=tev1_trials_orig(:,6)==-1 | (type_orig==1 & result_orig~=1)';
        %fixbreak_orig=fixbreak_orig|error_orig;
        fixbreak_orig=result_orig'==2|result_orig'==3|result_orig'==4;
        %fixbreak_orig=false(size(fixbreak_orig));
        
        sub(s).block(blocknum).bpseudo(1).result=result_orig;
        sub(s).block(blocknum).bpseudo(1).fixbreak=fixbreak_orig;
        
        sub(s).block(blocknum).bpseudo(2).result=result_orig;
        sub(s).block(blocknum).bpseudo(2).fixbreak=fixbreak_orig;
        
        for lead=1:length(sub(s).block(blocknum).b)
            lead
            for ev=1:length(event_ord)
                ev
                j=event_ord(ev);
                suj=suevent_ord(ev);
                for i=1:size(sub(s).block(blocknum).b(lead).su,2)
                    i
                    
                    trialtimes=trialtimes_orig(sub(s).block(blocknum).b(lead).su(suj,i).trials,:);
                    data=sub(s).block(blocknum).b(lead).su(suj,i).data;
                    event_happened=event_happened_orig(sub(s).block(blocknum).b(lead).su(suj,i).trials,:);
                    fixbreak=fixbreak_orig(sub(s).block(blocknum).b(lead).su(suj,i).trials);
                    tev1_trials=tev1_trials_orig(sub(s).block(blocknum).b(lead).su(suj,i).trials,:);
                    type=type_orig(sub(s).block(blocknum).b(lead).su(suj,i).trials);
                    pre=NaN*ones(length(type),1);
                    post=NaN*ones(length(type),1);
                    
                    switch j
                        case 1
                            pre(1)=-5;
                            pre(2:end)=trialtimes(1:end-1,2)-trialtimes(2:end,1)+preadj(ev);
                            
                            [row,col]=find(event_happened(:,ev+1:end));
                            nextev=accumarray(row,col,[size(event_happened,1),1],@min,NaN);
                            nextev=nextev+ev;
                            post=tev1_trials(sub2ind(size(tev1_trials),1:length(type),nextev'))'+postadj(ev)-tev1_trials(:,ev);
                            
                        case 8
                            [row,col]=find(event_happened(:,1:ev-1));
                            priorev=accumarray(row,col,[size(event_happened,1),1],@max,NaN);
                            pre=tev1_trials(sub2ind(size(tev1_trials),1:length(type),priorev'))'+preadj(ev)-tev1_trials(:,ev);
                            
                            post(1:end-1)=trialtimes(2:end,1)-trialtimes(1:end-1,2)+postadj(ev);
                            post(end)=5;
                        otherwise
                            [row,col]=find(event_happened(:,1:ev-1));
                            priorev=accumarray(row,col,[size(event_happened,1),1],@max,NaN);
                            pre=tev1_trials(sub2ind(size(tev1_trials),1:length(type),priorev'))'+preadj(ev)-tev1_trials(:,ev);
                            
                            [row,col]=find(event_happened(:,ev+1:end));
                            nextev=accumarray(row,col,[size(event_happened,1),1],@min,NaN);
                            nextev=nextev+ev;
                            post=tev1_trials(sub2ind(size(tev1_trials),1:length(type),nextev'))'+postadj(ev)-tev1_trials(:,ev);
                            
                    end
                    
                    
                    switch j
                        case 4
                            for k=[5]
                                if nnz(type==k)~=0
                                    [sub(s).block(blocknum).b(lead).su(suj,i).fr{k},sub(s).block(blocknum).b(lead).su(suj,i).fr_times{k},...
                                        sub(s).block(blocknum).b(lead).su(suj,i).fr_se{k}]=...
                                        psth_var(data(type'==k & ~fixbreak),[pre(~fixbreak & type'==k),post(~fixbreak & type'==k)],0.05,'n',[-2,2],2,-2:0.01:2);
                                end
                            end
                        case 9
                            [R,t,E,pmat,RR_ind]=psth_var_sig(data(~fixbreak),[pre(~fixbreak),post(~fixbreak)],0.05,'n',[-2,2],2,-2:0.01:2,type(~fixbreak),2,0.5,[1,2]);
                            sub(s).block(blocknum).b(lead).su(suj,i).pmat=pmat;
                            for k=2:3
                                if nnz(type==types(k))~=0
                                    sub(s).block(blocknum).b(lead).su(suj,i).fr{types(k)}=R(k-1,:);
                                    sub(s).block(blocknum).b(lead).su(suj,i).fr_times{types(k)}=t;
                                    sub(s).block(blocknum).b(lead).su(suj,i).fr_se{types(k)}=E(k-1,:);
                                    sub(s).block(blocknum).b(lead).su(suj,i).fr_ind{types(k)}=RR_ind{k-1};
                                end
                            end
                        otherwise
                            %if length(unique(type(~fixbreak))~=3)
                            %    continue
                            [R,t,E,pmat,RR_ind]=psth_var_sig(data(~fixbreak),[pre(~fixbreak),post(~fixbreak)],0.05,'n',[-2,2],2,-2:0.01:2,type(~fixbreak),2,0.5,[2,3]);
                            sub(s).block(blocknum).b(lead).su(suj,i).pmat=pmat;
                            for k=1:3
                                if nnz(type==types(k))~=0
                                    sub(s).block(blocknum).b(lead).su(suj,i).fr{types(k)}=R(k,:);
                                    sub(s).block(blocknum).b(lead).su(suj,i).fr_times{types(k)}=t;
                                    sub(s).block(blocknum).b(lead).su(suj,i).fr_se{types(k)}=E(k,:);
                                    sub(s).block(blocknum).b(lead).su(suj,i).fr_ind{types(k)}=RR_ind{k};
                                end
                            end
                            %end
                    end
                    presplit{s,blocknum,lead,ev,i}=pre;
                    postsplit{s,blocknum,lead,ev,i}=post;
                    fixbreakall{s,blocknum,lead,ev,i}=fixbreak;
                    typeall{s,blocknum,lead,ev,i}=type;
                    
                    preall{s,ev}=[preall{s,ev};pre(logical(event_happened(:,ev)))];
                    postall{s,ev}=[postall{s,ev};post(logical(event_happened(:,ev)))];
                    typeallall{s,ev}=[typeallall{s,ev};type(logical(event_happened(:,ev)))'];
                    fixbreakallall{s,ev}=[fixbreakallall{s,ev};fixbreak(logical(event_happened(:,ev)))];
                end
                
                toc
            end
            
        end
    end
end

%check that the same trials are used in all epochs for a cell, identification of
%cells where this is not true will be displayed


for s=[1,2,4,7,9,11]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            for i=1:size(sub(s).block(blocknum).b(lead).su,2)
                if nnz(squeeze(cellfun(@(x) nnz(~x), fixbreakall(s,blocknum,lead,[1,2,3,4,7],i)))~=nnz(~fixbreakall{s,blocknum,lead,1,i})*ones(5,1))>0
                    [s,blocknum,lead,i]
                end
                if nnz(squeeze(cellfun(@(x) nnz(~x), fixbreakall(s,blocknum,lead,[5,6],i)))~=nnz(~fixbreakall{s,blocknum,lead,5,i})*ones(2,1))>0
                    [s,blocknum,lead,i]
                end
            end
        end
    end
end

%% calculate PSTHs separated by trial type (requires chronux toolbox), pseudo trials
event_ord=[1,2,3,10,4,9,8];
suevent_ord=[1,2,3,9,4,8,7];
types=[1,4,5];

preadj=[0.5,0.05,0.05,0.05,0.05,0.05,0.05];
postadj=[0,0,-0.1,0,-0.1,0,0];

preall=cell(11,length(event_ord));
postall=preall;
typeallall=preall;
fixbreakallall=preall;


presplit=cell(11,5,2,length(event_ord),3);
postsplit=presplit;
typeall=presplit;
fixbreakall=presplit;tic
for s=[1,2,4,7,8,9,11]
    s
    for blocknum=1:length(sub(s).block)
        blocknum
        result_orig=sub(s).block(blocknum).bpseudo(1).result;
        type_orig=sub(s).block(blocknum).bpseudo(1).type;
        trialtimes_orig=sub(s).block(blocknum).bpseudo(1).t;
        tev1_trials_orig=sub(s).block(blocknum).bpseudo(1).tev1_trials;
        event_happened_orig=sub(s).block(blocknum).bpseudo(1).event_happened;
        event_happened_orig=[event_happened_orig,tev1_trials_orig(:,9)~=-1 & ~isnan(tev1_trials_orig(:,9)),tev1_trials_orig(:,10)~=-1 & ~isnan(tev1_trials_orig(:,10))];
        event_happened_orig(isnan(tev1_trials_orig))=0;
        event_happened_orig=event_happened_orig(:,event_ord);
        tev1_trials_orig=tev1_trials_orig(:,event_ord);
        %fixbreak_orig=tev1_trials_orig(:,4)<tev1_trials_orig(:,3) | isnan(tev1_trials_orig(:,2)) | isnan(tev1_trials_orig(:,3));
        %error_orig=tev1_trials_orig(:,6)==-1 & (type_orig'==1 & result_orig'~=1);
        %error_orig=result_orig'~=1;
        fixbreak_orig=sub(s).block(blocknum).bpseudo(1).fixbreak;
        
        
        for lead=1:length(sub(s).block(blocknum).b)
            lead
            for ev=1:length(event_ord)
                ev
                j=event_ord(ev);
                suj=suevent_ord(ev);
                for i=1:size(sub(s).block(blocknum).b(lead).su,2)
                    i
                    
                    trialtimes=trialtimes_orig(sub(s).block(blocknum).bpseudo(lead).su(suj,i).trials,:);
                    data=sub(s).block(blocknum).bpseudo(lead).su(suj,i).data;
                    event_happened=event_happened_orig(sub(s).block(blocknum).bpseudo(lead).su(suj,i).trials,:);
                    fixbreak=fixbreak_orig(sub(s).block(blocknum).bpseudo(lead).su(suj,i).trials);
                    tev1_trials=tev1_trials_orig(sub(s).block(blocknum).bpseudo(lead).su(suj,i).trials,:);
                    %error=error_orig(sub(s).block(blocknum).bpseudo(lead).su(suj,i).trials);
                    %ixbreak=fixbreak|error;
                    type=type_orig(sub(s).block(blocknum).bpseudo(lead).su(suj,i).trials);
                    pre=NaN*ones(length(type),1);
                    post=NaN*ones(length(type),1);
                    
                    switch j
                        
                        case 1
                            pre(1)=-5;
                            pre(2:end)=trialtimes(1:end-1,2)-trialtimes(2:end,1)+preadj(ev);
                            
                            [row,col]=find(event_happened(:,ev+1:end));
                            nextev=accumarray(row,col,[size(event_happened,1),1],@min,NaN);
                            nextev=nextev+ev;
                            post=tev1_trials(sub2ind(size(tev1_trials),1:length(type),nextev'))'+postadj(ev)-tev1_trials(:,ev);
                            
                        case 8
                            [row,col]=find(event_happened(:,1:ev-1));
                            priorev=accumarray(row,col,[size(event_happened,1),1],@max,NaN);
                            pre=tev1_trials(sub2ind(size(tev1_trials),1:length(type),priorev'))'+preadj(ev)-tev1_trials(:,ev);
                            
                            post(1:end-1)=trialtimes(2:end,1)-trialtimes(1:end-1,2)+postadj(ev);
                            post(end)=5;
                        otherwise
                            [row,col]=find(event_happened(:,1:ev-1));
                            priorev=accumarray(row,col,[size(event_happened,1),1],@max,NaN);
                            pre=tev1_trials(sub2ind(size(tev1_trials),1:length(type),priorev'))'+preadj(ev)-tev1_trials(:,ev);
                            
                            [row,col]=find(event_happened(:,ev+1:end));
                            nextev=accumarray(row,col,[size(event_happened,1),1],@min,NaN);
                            nextev=nextev+ev;
                            post=tev1_trials(sub2ind(size(tev1_trials),1:length(type),nextev'))'+postadj(ev)-tev1_trials(:,ev);
                            
                    end
                    
                    
                    
                    [R,t,E,pmat,RR_ind]=psth_var_sig(data(~fixbreak),[pre(~fixbreak),post(~fixbreak)],0.05,'n',[-2,2],2,-2:0.01:2,type(~fixbreak),2,0.5,[2,3]);
                    sub(s).block(blocknum).bpseudo(lead).su(suj,i).pmat=pmat;
                    for k=1:3
                        if nnz(type==types(k))~=0
                            sub(s).block(blocknum).bpseudo(lead).su(suj,i).fr{types(k)}=R(k,:);
                            sub(s).block(blocknum).bpseudo(lead).su(suj,i).fr_times{types(k)}=t;
                            sub(s).block(blocknum).bpseudo(lead).su(suj,i).fr_se{types(k)}=E(k,:);
                            sub(s).block(blocknum).bpseudo(lead).su(suj,i).fr_ind{types(k)}=RR_ind{k};
                        end
                    end
                    presplit{s,blocknum,lead,ev,i}=pre;
                    postsplit{s,blocknum,lead,ev,i}=post;
                    fixbreakall{s,blocknum,lead,ev,i}=fixbreak;
                    typeall{s,blocknum,lead,ev,i}=type;
                    
                    preall{s,ev}=[preall{s,ev};pre(logical(event_happened(:,ev)))];
                    postall{s,ev}=[postall{s,ev};post(logical(event_happened(:,ev)))];
                    typeallall{s,ev}=[typeallall{s,ev};type(logical(event_happened(:,ev)))'];
                    fixbreakallall{s,ev}=[fixbreakallall{s,ev};fixbreak(logical(event_happened(:,ev)))];
                end
                toc
            end
            
        end
    end
end
%% plot rasters with PSTH aligned to particular event

figure
colors=get(gca,'colororder');
%events=[1,2,3,8];
events=[1,2,3,9,4,8,7];
%events=[3,9,8,7];
%eventmap=[3,4,6,5];
%event_names={'Fixation','Target','Go','Turn'};
event_names={'Fixation','Target','Go','Movement','Turn Signal','Turn','Response'};

preadj=[0.5,0.1,0.1,0.1,0.1,0.1,0.1];
postadj=[0,0,-0.15,0,-0.15,0,0];

for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            if ~(s==1 && blocknum==1 && lead==1)
                figure
            end
            
            %            go=find(sub(s).block(blocknum).b(lead).type==1);
            %            turn=find(sub(s).block(blocknum).b(lead).type==4);
            %            swerve=find(sub(s).block(blocknum).b(lead).type==5);
            
            
            for i=1:size(sub(s).block(blocknum).b(lead).su,2)
                maxy=[];
                miny=[];
                for j=1:length(events)
                    if length(sub(s).block(blocknum).b(lead).su(events(j),i).strials)<100 % || ~isfield(sub(s).block(blocknum).b(lead).su(events(j),i),'fr')
                        continue
                    end
                    
                    event=events(j);
                    
                    result=sub(s).block(blocknum).b(1).result;
                    type=sub(s).block(blocknum).b(1).type;
                    tev1_trials=sub(s).block(blocknum).b(lead).tev1_trials;
                    tev1_trials(sub(s).block(blocknum).b(lead).type==4,4)=NaN;
                    tev1_trials=tev1_trials(:,event_ord);
                    %fixbreak=tev1_trials(:,4)<tev1_trials(:,3) | isnan(tev1_trials(:,2)) | isnan(tev1_trials(:,3));
                    %error=tev1_trials(:,6)==-1;
                    result(result==2 & type==4 & ~(tev1_trials(:,6)==-1 | isnan(tev1_trials(:,6)))')=1; %turn trials where a turn was identified but were marked incorrect
                    result(result==2 & type==5 & ~isnan(tev1_trials(:,5))' & ~(tev1_trials(:,6)==-1 | isnan(tev1_trials(:,6)))')=1; %swerve trials where a swerve signal was present and a turn was identified but were marked incorrect
                    result(result==1 & (type==4 | type==5) & tev1_trials(:,6)'==-1)=2; %turn or swerve trials marked correct but no turn identified
                    result(~isnan(tev1_trials(:,6)) & tev1_trials(:,6)~=-1 & tev1_trials(:,6)<tev1_trials(:,3))=4; %trials where a turn was identified before the go cue (i.e. incorrectly identified), ideally would go back and fix turn identification here
                    tev1_trials(result'~=3 & tev1_trials(:,4)<tev1_trials(:,3),4)=tev1_trials(result'~=3 & tev1_trials(:,4)<tev1_trials(:,3),3)+0.05;
                    fixbreak=result'==2|result'==3|result'==4;
                    
                    %fixbreak=fixbreak|error;
                    %fixbreak=false(size(fixbreak));
                    trials=zeros(1,length(sub(s).block(blocknum).b(lead).type));
                    trials(sub(s).block(blocknum).b(lead).su(events(j),i).trials)=true;
                    
                    go=find(sub(s).block(blocknum).b(lead).type'==1 & trials' & ~fixbreak);
                    turn=find(sub(s).block(blocknum).b(lead).type'==4 & trials' & ~fixbreak);
                    swerve=find(sub(s).block(blocknum).b(lead).type'==5 & trials' & ~fixbreak);
                    
                    %subplot(2*length(events),size(sub(s).block(blocknum).b(lead).su,2),2*(j-1)*size(sub(s).block(blocknum).b(lead).su,2)+i)
                    subplot(2*size(sub(s).block(blocknum).b(lead).su,2),length(events),2*(i-1)*length(events)+j)
                    hold on
                    switch event
                        case 4
                            swervelog=ismember(sub(s).block(blocknum).b(lead).su(events(j),i).strials,swerve);
                            swervetrials=sub(s).block(blocknum).b(lead).su(events(j),i).strials(swervelog);
                            [swervetrialsu,~,~]=unique(swervetrials);
                            swervespks=sub(s).block(blocknum).b(lead).su(events(j),i).stimes_align(swervelog);
                            swervetrialslog=typeall{s,blocknum,lead,j,i}'==5 & ~fixbreakall{s,blocknum,lead,j,i};
                            swervecounts=histcounts(swervetrials,[swervetrialsu-0.5;max(swervetrialsu)+0.5]);
                            swervecounts=mat2cell(swervecounts',ones(1,nnz(swervetrialslog)));
                            preswerve=mat2cell(presplit{s,blocknum,lead,j,i}(swervetrialslog),ones(1,nnz(swervetrialslog)));
                            preswerve=cellfun(@(x,y) repmat(x,1,y)',preswerve,swervecounts,'uniformoutput',0);
                            postswerve=mat2cell(postsplit{s,blocknum,lead,j,i}(swervetrialslog),ones(1,nnz(swervetrialslog)));
                            postswerve=cellfun(@(x,y) repmat(x,1,y)',postswerve,swervecounts,'uniformoutput',0);
                            swervespks(swervespks<vertcat(preswerve{:}))=NaN;
                            swervespks(swervespks>vertcat(postswerve{:}))=NaN;
                            [~,~,swervetrials]=unique(swervetrials);
                            
                            plot(swervespks*1000,swervetrials,'.','MarkerSize',3,'MarkerEdgeColor',colors(3,:)) %plot rasters
                        case 8
                            turnlog=ismember(sub(s).block(blocknum).b(lead).su(events(j),i).strials,turn);
                            turntrials=sub(s).block(blocknum).b(lead).su(events(j),i).strials(turnlog);
                            [turntrialsu,~,~]=unique(turntrials);
                            turnspks=sub(s).block(blocknum).b(lead).su(events(j),i).stimes_align(turnlog);
                            turntrialslog=typeall{s,blocknum,lead,j,i}'==4 & ~fixbreakall{s,blocknum,lead,j,i};
                            turncounts=histcounts(turntrials,[turntrialsu-0.5;max(turntrialsu)+0.5]);
                            turncounts=mat2cell(turncounts',ones(1,nnz(turntrialslog)));
                            preturn=mat2cell(presplit{s,blocknum,lead,j,i}(turntrialslog),ones(1,nnz(turntrialslog)));
                            preturn=cellfun(@(x,y) repmat(x,1,y)',preturn,turncounts,'uniformoutput',0);
                            postturn=mat2cell(postsplit{s,blocknum,lead,j,i}(turntrialslog),ones(1,nnz(turntrialslog)));
                            postturn=cellfun(@(x,y) repmat(x,1,y)',postturn,turncounts,'uniformoutput',0);
                            turnspks(turnspks<vertcat(preturn{:}))=NaN;
                            turnspks(turnspks>vertcat(postturn{:}))=NaN;
                            [~,~,turntrials]=unique(turntrials);
                            
                            swervelog=ismember(sub(s).block(blocknum).b(lead).su(events(j),i).strials,swerve);
                            swervetrials=sub(s).block(blocknum).b(lead).su(events(j),i).strials(swervelog);
                            [swervetrialsu,~,~]=unique(swervetrials);
                            swervespks=sub(s).block(blocknum).b(lead).su(events(j),i).stimes_align(swervelog);
                            swervetrialslog=typeall{s,blocknum,lead,j,i}'==5 & ~fixbreakall{s,blocknum,lead,j,i};
                            swervecounts=histcounts(swervetrials,[swervetrialsu-0.5;max(swervetrialsu)+0.5]);
                            swervecounts=mat2cell(swervecounts',ones(1,nnz(swervetrialslog)));
                            preswerve=mat2cell(presplit{s,blocknum,lead,j,i}(swervetrialslog),ones(1,nnz(swervetrialslog)));
                            preswerve=cellfun(@(x,y) repmat(x,1,y)',preswerve,swervecounts,'uniformoutput',0);
                            postswerve=mat2cell(postsplit{s,blocknum,lead,j,i}(swervetrialslog),ones(1,nnz(swervetrialslog)));
                            postswerve=cellfun(@(x,y) repmat(x,1,y)',postswerve,swervecounts,'uniformoutput',0);
                            swervespks(swervespks<vertcat(preswerve{:}))=NaN;
                            swervespks(swervespks>vertcat(postswerve{:}))=NaN;
                            [~,~,swervetrials]=unique(swervetrials);
                            swervetrials=swervetrials+turntrials(end);
                            
                            plot(turnspks*1000,turntrials,'.','MarkerSize',3,'MarkerEdgeColor',colors(2,:)) %plot rasters
                            plot(swervespks*1000,swervetrials,'.','MarkerSize',3,'MarkerEdgeColor',colors(3,:)) %plot rasters
                        otherwise
                            golog=ismember(sub(s).block(blocknum).b(lead).su(events(j),i).strials,go);
                            gotrials=sub(s).block(blocknum).b(lead).su(events(j),i).strials(golog);
                            [gotrialsu,~,~]=unique(gotrials);
                            gospks=sub(s).block(blocknum).b(lead).su(events(j),i).stimes_align(golog);
                            gotrialslog=typeall{s,blocknum,lead,j,i}'==1 & ~fixbreakall{s,blocknum,lead,j,i};
                            gocounts=histcounts(gotrials,[gotrialsu-0.5;max(gotrialsu)+0.5]);
                            gocounts=mat2cell(gocounts',ones(1,nnz(gotrialslog)));
                            prego=mat2cell(presplit{s,blocknum,lead,j,i}(gotrialslog),ones(1,nnz(gotrialslog)));
                            prego=cellfun(@(x,y) repmat(x,1,y)',prego,gocounts,'uniformoutput',0);
                            postgo=mat2cell(postsplit{s,blocknum,lead,j,i}(gotrialslog),ones(1,nnz(gotrialslog)));
                            postgo=cellfun(@(x,y) repmat(x,1,y)',postgo,gocounts,'uniformoutput',0);
                            gospks(gospks<vertcat(prego{:}))=NaN;
                            gospks(gospks>vertcat(postgo{:}))=NaN;
                            [~,~,gotrials]=unique(gotrials);
                            
                            turnlog=ismember(sub(s).block(blocknum).b(lead).su(events(j),i).strials,turn);
                            turntrials=sub(s).block(blocknum).b(lead).su(events(j),i).strials(turnlog);
                            [turntrialsu,~,~]=unique(turntrials);
                            turnspks=sub(s).block(blocknum).b(lead).su(events(j),i).stimes_align(turnlog);
                            turntrialslog=typeall{s,blocknum,lead,j,i}'==4 & ~fixbreakall{s,blocknum,lead,j,i};
                            turncounts=histcounts(turntrials,[turntrialsu-0.5;max(turntrialsu)+0.5]);
                            turncounts=mat2cell(turncounts',ones(1,nnz(turntrialslog)));
                            preturn=mat2cell(presplit{s,blocknum,lead,j,i}(turntrialslog),ones(1,nnz(turntrialslog)));
                            preturn=cellfun(@(x,y) repmat(x,1,y)',preturn,turncounts,'uniformoutput',0);
                            postturn=mat2cell(postsplit{s,blocknum,lead,j,i}(turntrialslog),ones(1,nnz(turntrialslog)));
                            postturn=cellfun(@(x,y) repmat(x,1,y)',postturn,turncounts,'uniformoutput',0);
                            turnspks(turnspks<vertcat(preturn{:}))=NaN;
                            turnspks(turnspks>vertcat(postturn{:}))=NaN;
                            [~,~,turntrials]=unique(turntrials);
                            turntrials=turntrials+gotrials(end);
                            
                            swervelog=ismember(sub(s).block(blocknum).b(lead).su(events(j),i).strials,swerve);
                            swervetrials=sub(s).block(blocknum).b(lead).su(events(j),i).strials(swervelog);
                            [swervetrialsu,~,~]=unique(swervetrials);
                            swervespks=sub(s).block(blocknum).b(lead).su(events(j),i).stimes_align(swervelog);
                            swervetrialslog=typeall{s,blocknum,lead,j,i}'==5 & ~fixbreakall{s,blocknum,lead,j,i};
                            swervecounts=histcounts(swervetrials,[swervetrialsu-0.5;max(swervetrialsu)+0.5]);
                            swervecounts=mat2cell(swervecounts',ones(1,nnz(swervetrialslog)));
                            preswerve=mat2cell(presplit{s,blocknum,lead,j,i}(swervetrialslog),ones(1,nnz(swervetrialslog)));
                            preswerve=cellfun(@(x,y) repmat(x,1,y)',preswerve,swervecounts,'uniformoutput',0);
                            postswerve=mat2cell(postsplit{s,blocknum,lead,j,i}(swervetrialslog),ones(1,nnz(swervetrialslog)));
                            postswerve=cellfun(@(x,y) repmat(x,1,y)',postswerve,swervecounts,'uniformoutput',0);
                            swervespks(swervespks<vertcat(preswerve{:}))=NaN;
                            swervespks(swervespks>vertcat(postswerve{:}))=NaN;
                            [~,~,swervetrials]=unique(swervetrials);
                            swervetrials=swervetrials+turntrials(end);
                            
                            plot(gospks*1000,gotrials,'.','MarkerSize',3,'MarkerEdgeColor',colors(1,:)) %plot rasters
                            plot(turnspks*1000, turntrials,'.','MarkerSize',3,'MarkerEdgeColor',colors(2,:)) %plot rasters
                            plot(swervespks*1000,swervetrials,'.','MarkerSize',3,'MarkerEdgeColor',colors(3,:)) %plot rasters
                            
                            
                    end
                    
                    %plot(sub(s).block(blocknum).b(lead).su(events(j),i).stimes_align*1000,sub(s).block(blocknum).b(lead).su(events(j),i).strials,'.','MarkerSize',1,'MarkerEdgeColor','k') %plot rasters
                    %xlim([-2000,2000])
                    
                    if event~=9
                        %xlim(1000*[median(preall{s,j})+preadj(j),median(postall{s,j})+postadj(j)])
                        xlim(1000*[median(preall{s,j}),median(postall{s,j})])
                    else
                        swervelog2=[];
                        type2=[];
                        fixbreak=[];
                        for blocknum2=1:length(sub(s).block)
                            for lead2=1:length(sub(s).block(blocknum2).b)
                                for i2=1:size(sub(s).block(blocknum2).b(lead2).su,2)
                                    tev1_trials=sub(s).block(blocknum2).b(lead2).tev1_trials(:,10);
                                    tev1_trials=tev1_trials(sub(s).block(blocknum2).b(lead2).su(9,i2).trials,:);
                                    type2=[type2;typeall{s,blocknum2,lead2,j,i2}'];
                                    fixbreak=[fixbreak;fixbreakall{s,blocknum2,lead2,j,i2}];
                                end
                            end
                        end
                        postsort=sort(postall{s,j}(type2==5 & ~fixbreak))+postadj(j);
                        %xlim(1000*[median(preall{s,j})+preadj(j),postsort(round(length(postsort)/4))])
                        xlim(1000*[median(preall{s,j}),postsort(round(length(postsort)/4))])
                    end
                    ylim([0,swervetrials(end)+1])
                    %ylim([0,max(sub(s).block(blocknum).b(lead).su(events(j),i).strials)+1])
                    %qual=sub(s).block(blocknum).b(lead).quality(i);
                    qual=0;
                    if j==1
                        %qual=sub(s).block(blocknum).b(lead).quality(i);
                        qual=0;
                        if qual==0
                            qual='hash';
                        elseif qual==5
                            qual='multi';
                        else
                            qual=num2str(qual);
                        end
                        
                        text(-0.7,-0.5,['Unit ',num2str(i),' Quality: ',qual],'units','normalize','rotation',90,'horizontalalignment','center')
                        ylabel('Trial')
                    else
                        set(gca,'YTickLabel',{})
                    end
                    
                    if i==1
                        title(event_names{j})
                    end
                    set(gca,'box','off','XTickLabel',{},'FontSize',12);
                    %subplot(2*length(events),size(sub(s).block(blocknum).b(lead).su,2),2*(j-1)*size(sub(s).block(blocknum).b(lead).su,2)+size(sub(s).block(blocknum).b(lead).su,2)+i)
                    frplot(j)=subplot(2*size(sub(s).block(blocknum).b(lead).su,2),length(events),2*(i-1)*length(events)+length(events)+j);
                    patchymat=[];
                    %if j==4
                    if j==6
                        
                        col=2;
                        for k=[4,5]
                            patchx=[sub(s).block(blocknum).b(lead).su(events(j),i).fr_times{k}*1000,fliplr(sub(s).block(blocknum).b(lead).su(events(j),i).fr_times{k}*1000)];
                            patchy=[sub(s).block(blocknum).b(lead).su(events(j),i).fr{k}-sub(s).block(blocknum).b(lead).su(events(j),i).fr_se{k},fliplr(sub(s).block(blocknum).b(lead).su(events(j),i).fr{k}+sub(s).block(blocknum).b(lead).su(events(j),i).fr_se{k})];
                            err(k)=patch(patchx(~isnan(patchy)),patchy(~isnan(patchy)),mean([colors(col,:);1,1,1]),'EdgeColor','none','FaceAlpha',0.5);
                            
                            hold on
                            
                            l(k)=plot(sub(s).block(blocknum).b(lead).su(events(j),i).fr_times{k}*1000,sub(s).block(blocknum).b(lead).su(events(j),i).fr{k},'Color',colors(col,:));
                            col=col+1;
                            patchymat=[patchymat;patchy];
                        end
                    elseif j==5
                        col=3;
                        for k=[5]
                            patchx=[sub(s).block(blocknum).b(lead).su(events(j),i).fr_times{k}*1000,fliplr(sub(s).block(blocknum).b(lead).su(events(j),i).fr_times{k}*1000)];
                            patchy=[sub(s).block(blocknum).b(lead).su(events(j),i).fr{k}-sub(s).block(blocknum).b(lead).su(events(j),i).fr_se{k},fliplr(sub(s).block(blocknum).b(lead).su(events(j),i).fr{k}+sub(s).block(blocknum).b(lead).su(events(j),i).fr_se{k})];
                            err(k)=patch(patchx(~isnan(patchy)),patchy(~isnan(patchy)),mean([colors(col,:);1,1,1]),'EdgeColor','none','FaceAlpha',0.5);
                            
                            hold on
                            
                            l(k)=plot(sub(s).block(blocknum).b(lead).su(events(j),i).fr_times{k}*1000,sub(s).block(blocknum).b(lead).su(events(j),i).fr{k},'Color',colors(col,:));
                            col=col+1;
                            patchymat=[patchymat;patchy];
                        end
                    else
                        col=1;
                        for k=[1,4,5]
                            patchx=[sub(s).block(blocknum).b(lead).su(events(j),i).fr_times{k}*1000,fliplr(sub(s).block(blocknum).b(lead).su(events(j),i).fr_times{k}*1000)];
                            patchy=[sub(s).block(blocknum).b(lead).su(events(j),i).fr{k}-sub(s).block(blocknum).b(lead).su(events(j),i).fr_se{k},fliplr(sub(s).block(blocknum).b(lead).su(events(j),i).fr{k}+sub(s).block(blocknum).b(lead).su(events(j),i).fr_se{k})];
                            err(k)=patch(patchx(~isnan(patchy)),patchy(~isnan(patchy)),mean([colors(col,:);1,1,1]),'EdgeColor','none','FaceAlpha',0.5);
                            
                            hold on
                            
                            l(k)=plot(sub(s).block(blocknum).b(lead).su(events(j),i).fr_times{k}*1000,sub(s).block(blocknum).b(lead).su(events(j),i).fr{k},'Color',colors(col,:));
                            col=col+1;
                            patchymat=[patchymat;patchy];
                        end
                    end
                    %xlim([-2000,2000])
                    ylims=get(gca,'YLim');
                    if event~=9
                        %xlim(1000*[median(preall{s,j})+preadj(j),min([median(postall{s,j})+postadj(j),2])])
                        xlim(1000*[median(preall{s,j}),min([median(postall{s,j}),2])])
                        %boundinds=discretize([median(preall{s,j})+preadj(j),min([median(postall{s,j})+postadj(j),2])],patchx(1:(length(patchx)/2))./1000);
                        boundinds=discretize([median(preall{s,j}),min([median(postall{s,j}),2])],patchx(1:(length(patchx)/2))./1000);
                        maxy=max(max(patchymat(:,(length(patchymat)-boundinds(2)):(length(patchymat)-boundinds(1))),[],2));
                        miny=min(min(patchymat(:,boundinds(1):boundinds(2)),[],2));
                    else
                        swervelog2=[];
                        type2=[];
                        fixbreak=[];
                        for blocknum2=1:length(sub(s).block)
                            for lead2=1:length(sub(s).block(blocknum2).b)
                                for i2=1:size(sub(s).block(blocknum2).b(lead2).su,2)
                                    tev1_trials=sub(s).block(blocknum2).b(lead2).tev1_trials(:,10);
                                    tev1_trials=tev1_trials(sub(s).block(blocknum2).b(lead2).su(9,i2).trials,:);
                                    type2=[type2;typeall{s,blocknum2,lead2,j,i2}'];
                                    fixbreak=[fixbreak;fixbreakall{s,blocknum2,lead2,j,i2}];
                                end
                            end
                        end
                        postsort=sort(postall{s,j}(type2==5 & ~fixbreak))+postadj(j);
                        %xlim(1000*[median(preall{s,j})+preadj(j),min([postsort(round(length(postsort)/4)),2])])
                        %boundinds=discretize([median(preall{s,j})+preadj(j),min([postsort(round(length(postsort)/4)),2])],patchx(1:(length(patchx)/2))./1000);
                        xlim(1000*[median(preall{s,j}),min([postsort(round(length(postsort)/4)),2])])
                        boundinds=discretize([median(preall{s,j}),min([postsort(round(length(postsort)/4)),2])],patchx(1:(length(patchx)/2))./1000);
                        maxy=max(max(patchymat(:,(length(patchymat)-boundinds(2)):(length(patchymat)-boundinds(1))),[],2));
                        miny=min(min(patchymat(:,boundinds(1):boundinds(2)),[],2));
                    end
                    %if event~=4
                    %    sig=(sub(s).block(blocknum).b(lead).su(events(j),i).pmat<0.0001);
                    %    sig=double(sig);
                    %    sig(sig==0)=NaN;
                    %    plot(sub(s).block(blocknum).b(lead).su(events(j),i).fr_times{k}*1000,sig*ylims(2)*0.9,'LineWidth',2,'color','k');
                    %end
                    %ylims(1)=0;
                    %ylims(2)=1.1*max(ymax);
                    %if j==1
                    %    ymax=ylims(2);
                    %else
                    %    ylims(2)=ymax;
                    %end
                    %set(gca,'YLim',ylims,'FontSize',8);
                    %title(['Single Unit ',num2str(i)])
                    if i==size(sub(s).block(blocknum).b(lead).su,2)
                        xlabel('Time (ms)')
                        %if j==length(events)
                        %    leg=legend(err([1,4,5]),'Go','Turn','Swerve');
                        %    set(leg,'box','off','Position',[0.47,-0.01,0.1,0.05],'units','normalized','orientation','horizontal')
                        %    if size(sub(s).block(blocknum).b(lead).su,2)>1
                        %        %tightfig;
                        %    end
                        %end
                    else
                        %set(gca,'XTickLabel',{})
                    end
                    if j==1
                        ylabel('FR (Hz)')
                    else
                        set(gca,'YTickLabel',{})
                    end
                    set(gca,'box','off','FontSize',12)
                    %if i~=1
                    %
                    %end
                end
                ylims(1)=0.8*min(miny);
                ylims(2)=1.2*max(maxy);
                for j=1:length(events)
                    set(frplot(j),'YLim',ylims);
                end
            end
            suptitle(['Subject ',num2str(s),' Block ',num2str(blocknum),' Lead ',num2str(lead)])
            set(gcf,'Units','inches','Position',[3,3,12,3+2*size(sub(s).block(blocknum).b(lead).su,2)])
            print(['/Users/Dennis/Documents/MATLAB/Mogilner_Lab/interrupted_reach/Subject_',num2str(s),'_Block_',num2str(blocknum),'_Lead_',num2str(lead)],'-depsc')
        end
    end
end

%% plot rasters with PSTH aligned to particular event pseudo trials

figure
colors=get(gca,'colororder');
%events=[1,2,3,8];
events=[1,2,3,9,4,8,7];
%events=[3,9,8,7];
%eventmap=[3,4,6,5];
%event_names={'Fixation','Target','Go','Turn'};
event_names={'Fixation','Instruction','Go','Movement Onset','Turn Cue','Turn Onset','Feedback'};

preadj=[0.5,0.05,0.05,0.05,0.05,0.05,0.05];
postadj=[0,0,-0.1,0,-0.1,0,0];

for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            if ~(s==1 && blocknum==1 && lead==1)
                figure
            end
            
            %            go=find(sub(s).block(blocknum).bpseudo(lead).type==1);
            %            turn=find(sub(s).block(blocknum).bpseudo(lead).type==4);
            %            swerve=find(sub(s).block(blocknum).bpseudo(lead).type==5);
            
            
            for i=1:size(sub(s).block(blocknum).bpseudo(lead).su,2)
                maxy=[];
                miny=[];
                for j=1:length(events)
                    if length(sub(s).block(blocknum).bpseudo(lead).su(events(j),i).strials)<100 % || ~isfield(sub(s).block(blocknum).bpseudo(lead).su(events(j),i),'fr')
                        continue
                    end
                    
                    event=events(j);
                    
                    result=sub(s).block(blocknum).bpseudo(1).result;
                    type=sub(s).block(blocknum).bpseudo(1).type;
                    tev1_trials=sub(s).block(blocknum).bpseudo(lead).tev1_trials;
                    tev1_trials=tev1_trials(:,event_ord);
                    %fixbreak=tev1_trials(:,4)<tev1_trials(:,3) | isnan(tev1_trials(:,2)) | isnan(tev1_trials(:,3));
                    %error=tev1_trials(:,6)==-1;
                    fixbreak=sub(s).block(blocknum).bpseudo(1).fixbreak;
                    
                    %fixbreak=fixbreak|error;
                    %fixbreak=false(size(fixbreak));
                    trials=zeros(1,length(sub(s).block(blocknum).bpseudo(lead).type));
                    trials(sub(s).block(blocknum).bpseudo(lead).su(events(j),i).trials)=true;
                    
                    go=find(sub(s).block(blocknum).bpseudo(lead).type'==1 & trials' & ~fixbreak);
                    turn=find(sub(s).block(blocknum).bpseudo(lead).type'==4 & trials' & ~fixbreak);
                    swerve=find(sub(s).block(blocknum).bpseudo(lead).type'==5 & trials' & ~fixbreak);
                    
                    %subplot(2*length(events),size(sub(s).block(blocknum).bpseudo(lead).su,2),2*(j-1)*size(sub(s).block(blocknum).bpseudo(lead).su,2)+i)
                    subplot(2*size(sub(s).block(blocknum).bpseudo(lead).su,2),length(events),2*(i-1)*length(events)+j)
                    hold on
                    golog=ismember(sub(s).block(blocknum).bpseudo(lead).su(events(j),i).strials,go);
                    gotrials=sub(s).block(blocknum).bpseudo(lead).su(events(j),i).strials(golog);
                    [gotrialsu,~,~]=unique(gotrials);
                    gospks=sub(s).block(blocknum).bpseudo(lead).su(events(j),i).stimes_align(golog);
                    gotrialslog=typeall{s,blocknum,lead,j,i}'==1 & ~fixbreakall{s,blocknum,lead,j,i};
                    gocounts=histcounts(gotrials,[gotrialsu-0.5;max(gotrialsu)+0.5]);
                    gocounts=mat2cell(gocounts',ones(1,nnz(gotrialslog)));
                    prego=mat2cell(presplit{s,blocknum,lead,j,i}(gotrialslog),ones(1,nnz(gotrialslog)));
                    prego=cellfun(@(x,y) repmat(x,1,y)',prego,gocounts,'uniformoutput',0);
                    postgo=mat2cell(postsplit{s,blocknum,lead,j,i}(gotrialslog),ones(1,nnz(gotrialslog)));
                    postgo=cellfun(@(x,y) repmat(x,1,y)',postgo,gocounts,'uniformoutput',0);
                    gospks(gospks<vertcat(prego{:}))=NaN;
                    gospks(gospks>vertcat(postgo{:}))=NaN;
                    [~,~,gotrials]=unique(gotrials);
                    
                    turnlog=ismember(sub(s).block(blocknum).bpseudo(lead).su(events(j),i).strials,turn);
                    turntrials=sub(s).block(blocknum).bpseudo(lead).su(events(j),i).strials(turnlog);
                    [turntrialsu,~,~]=unique(turntrials);
                    turnspks=sub(s).block(blocknum).bpseudo(lead).su(events(j),i).stimes_align(turnlog);
                    turntrialslog=typeall{s,blocknum,lead,j,i}'==4 & ~fixbreakall{s,blocknum,lead,j,i};
                    turncounts=histcounts(turntrials,[turntrialsu-0.5;max(turntrialsu)+0.5]);
                    turncounts=mat2cell(turncounts',ones(1,nnz(turntrialslog)));
                    preturn=mat2cell(presplit{s,blocknum,lead,j,i}(turntrialslog),ones(1,nnz(turntrialslog)));
                    preturn=cellfun(@(x,y) repmat(x,1,y)',preturn,turncounts,'uniformoutput',0);
                    postturn=mat2cell(postsplit{s,blocknum,lead,j,i}(turntrialslog),ones(1,nnz(turntrialslog)));
                    postturn=cellfun(@(x,y) repmat(x,1,y)',postturn,turncounts,'uniformoutput',0);
                    turnspks(turnspks<vertcat(preturn{:}))=NaN;
                    turnspks(turnspks>vertcat(postturn{:}))=NaN;
                    [~,~,turntrials]=unique(turntrials);
                    turntrials=turntrials+gotrials(end);
                    
                    swervelog=ismember(sub(s).block(blocknum).bpseudo(lead).su(events(j),i).strials,swerve);
                    swervetrials=sub(s).block(blocknum).bpseudo(lead).su(events(j),i).strials(swervelog);
                    [swervetrialsu,~,~]=unique(swervetrials);
                    swervespks=sub(s).block(blocknum).bpseudo(lead).su(events(j),i).stimes_align(swervelog);
                    swervetrialslog=typeall{s,blocknum,lead,j,i}'==5 & ~fixbreakall{s,blocknum,lead,j,i};
                    swervecounts=histcounts(swervetrials,[swervetrialsu-0.5;max(swervetrialsu)+0.5]);
                    swervecounts=mat2cell(swervecounts',ones(1,nnz(swervetrialslog)));
                    preswerve=mat2cell(presplit{s,blocknum,lead,j,i}(swervetrialslog),ones(1,nnz(swervetrialslog)));
                    preswerve=cellfun(@(x,y) repmat(x,1,y)',preswerve,swervecounts,'uniformoutput',0);
                    postswerve=mat2cell(postsplit{s,blocknum,lead,j,i}(swervetrialslog),ones(1,nnz(swervetrialslog)));
                    postswerve=cellfun(@(x,y) repmat(x,1,y)',postswerve,swervecounts,'uniformoutput',0);
                    swervespks(swervespks<vertcat(preswerve{:}))=NaN;
                    swervespks(swervespks>vertcat(postswerve{:}))=NaN;
                    [~,~,swervetrials]=unique(swervetrials);
                    swervetrials=swervetrials+turntrials(end);
                    
                    plot(gospks*1000,gotrials,'.','MarkerSize',3,'MarkerEdgeColor',colors(1,:)) %plot rasters
                    plot(turnspks*1000, turntrials,'.','MarkerSize',3,'MarkerEdgeColor',colors(2,:)) %plot rasters
                    plot(swervespks*1000,swervetrials,'.','MarkerSize',3,'MarkerEdgeColor',colors(3,:)) %plot rasters
                    
                    
                    %plot(sub(s).block(blocknum).bpseudo(lead).su(events(j),i).stimes_align*1000,sub(s).block(blocknum).bpseudo(lead).su(events(j),i).strials,'.','MarkerSize',1,'MarkerEdgeColor','k') %plot rasters
                    %xlim([-2000,2000])
                    
                    %if event~=9
                    %xlim(1000*[median(preall{s,j})+preadj(j),median(postall{s,j})+postadj(j)])
                    xlim(1000*[median(preall{s,j}),median(postall{s,j})])
                    %else
                    %    swervelog2=[];
                    %    type2=[];
                    %    fixbreak=[];
                    %    for blocknum2=1:length(sub(s).block)
                    %        for lead2=1:length(sub(s).block(blocknum2).b)
                    %            for i2=1:size(sub(s).block(blocknum2).bpseudo(lead2).su,2)
                    %                tev1_trials=sub(s).block(blocknum2).bpseudo(lead2).tev1_trials(:,10);
                    %                tev1_trials=tev1_trials(sub(s).block(blocknum2).bpseudo(lead2).su(9,i2).trials,:);
                    %                type2=[type2;typeall{s,blocknum2,lead2,j,i2}'];
                    %                fixbreak=[fixbreak;fixbreakall{s,blocknum2,lead2,j,i2}];
                    %            end
                    %        end
                    %    end
                    %    postsort=sort(postall{s,j}(type2==5 & ~fixbreak))+postadj(j);
                    %xlim(1000*[median(preall{s,j})+preadj(j),postsort(round(length(postsort)/4))])
                    %    xlim(1000*[median(preall{s,j}),postsort(round(length(postsort)/4))])
                    %end
                    ylim([0,swervetrials(end)+1])
                    %ylim([0,max(sub(s).block(blocknum).bpseudo(lead).su(events(j),i).strials)+1])
                    %qual=sub(s).block(blocknum).bpseudo(lead).quality(i);
                    qual=0;
                    if j==1
                        qual=sub(s).block(blocknum).b(lead).quality(i);
                        %qual=0;
                        if qual==0
                            qual='hash';
                        elseif qual==5
                            qual='multi';
                        else
                            qual=num2str(qual);
                        end
                        
                        text(-0.7,-0.5,['Unit ',num2str(i),' Quality: ',qual],'units','normalize','rotation',90,'horizontalalignment','center')
                        ylabel('Trial')
                    else
                        set(gca,'YTickLabel',{})
                    end
                    
                    if i==1
                        title(event_names{j})
                    end
                    set(gca,'box','off','XTickLabel',{},'FontSize',12);
                    %subplot(2*length(events),size(sub(s).block(blocknum).bpseudo(lead).su,2),2*(j-1)*size(sub(s).block(blocknum).bpseudo(lead).su,2)+size(sub(s).block(blocknum).bpseudo(lead).su,2)+i)
                    frplot(j)=subplot(2*size(sub(s).block(blocknum).bpseudo(lead).su,2),length(events),2*(i-1)*length(events)+length(events)+j);
                    patchymat=[];
                    
                    col=1;
                    for k=[1,4,5]
                        patchx=[sub(s).block(blocknum).bpseudo(lead).su(events(j),i).fr_times{k}*1000,fliplr(sub(s).block(blocknum).bpseudo(lead).su(events(j),i).fr_times{k}*1000)];
                        patchy=[sub(s).block(blocknum).bpseudo(lead).su(events(j),i).fr{k}-sub(s).block(blocknum).bpseudo(lead).su(events(j),i).fr_se{k},fliplr(sub(s).block(blocknum).bpseudo(lead).su(events(j),i).fr{k}+sub(s).block(blocknum).bpseudo(lead).su(events(j),i).fr_se{k})];
                        err(k)=patch(patchx(~isnan(patchy)),patchy(~isnan(patchy)),mean([colors(col,:);1,1,1]),'EdgeColor','none','FaceAlpha',0.5);
                        
                        hold on
                        
                        l(k)=plot(sub(s).block(blocknum).bpseudo(lead).su(events(j),i).fr_times{k}*1000,sub(s).block(blocknum).bpseudo(lead).su(events(j),i).fr{k},'Color',colors(col,:));
                        col=col+1;
                        patchymat=[patchymat;patchy];
                    end
                    %xlim([-2000,2000])
                    ylims=get(gca,'YLim');
                    %if event~=9
                    %xlim(1000*[median(preall{s,j})+preadj(j),min([median(postall{s,j})+postadj(j),2])])
                    xlim(1000*[median(preall{s,j}),min([median(postall{s,j}),2])])
                    %boundinds=discretize([median(preall{s,j})+preadj(j),min([median(postall{s,j})+postadj(j),2])],patchx(1:(length(patchx)/2))./1000);
                    boundinds=discretize([median(preall{s,j}),min([median(postall{s,j}),2])],patchx(1:(length(patchx)/2))./1000);
                    maxy=max(max(patchymat(:,(length(patchymat)-boundinds(2)):(length(patchymat)-boundinds(1))),[],2));
                    miny=min(min(patchymat(:,boundinds(1):boundinds(2)),[],2));
                    %else
                    %    swervelog2=[];
                    %    type2=[];
                    %    fixbreak=[];
                    %    for blocknum2=1:length(sub(s).block)
                    %        for lead2=1:length(sub(s).block(blocknum2).b)
                    %            for i2=1:size(sub(s).block(blocknum2).bpseudo(lead2).su,2)
                    %                tev1_trials=sub(s).block(blocknum2).bpseudo(lead2).tev1_trials(:,10);
                    %                tev1_trials=tev1_trials(sub(s).block(blocknum2).bpseudo(lead2).su(9,i2).trials,:);
                    %                type2=[type2;typeall{s,blocknum2,lead2,j,i2}'];
                    %                fixbreak=[fixbreak;fixbreakall{s,blocknum2,lead2,j,i2}];
                    %            end
                    %        end
                    %    end
                    %    postsort=sort(postall{s,j}(type2==5 & ~fixbreak))+postadj(j);
                    %    %xlim(1000*[median(preall{s,j})+preadj(j),min([postsort(round(length(postsort)/4)),2])])
                    %    %boundinds=discretize([median(preall{s,j})+preadj(j),min([postsort(round(length(postsort)/4)),2])],patchx(1:(length(patchx)/2))./1000);
                    %    xlim(1000*[median(preall{s,j}),min([postsort(round(length(postsort)/4)),2])])
                    %    boundinds=discretize([median(preall{s,j}),min([postsort(round(length(postsort)/4)),2])],patchx(1:(length(patchx)/2))./1000);
                    %    maxy=max(max(patchymat(:,(length(patchymat)-boundinds(2)):(length(patchymat)-boundinds(1))),[],2));
                    %    miny=min(min(patchymat(:,boundinds(1):boundinds(2)),[],2));
                    %end
                    %if event~=4
                    %    sig=(sub(s).block(blocknum).bpseudo(lead).su(events(j),i).pmat<0.0001);
                    %    sig=double(sig);
                    %    sig(sig==0)=NaN;
                    %    plot(sub(s).block(blocknum).bpseudo(lead).su(events(j),i).fr_times{k}*1000,sig*ylims(2)*0.9,'LineWidth',2,'color','k');
                    %end
                    %ylims(1)=0;
                    %ylims(2)=1.1*max(ymax);
                    %if j==1
                    %    ymax=ylims(2);
                    %else
                    %    ylims(2)=ymax;
                    %end
                    %set(gca,'YLim',ylims,'FontSize',8);
                    %title(['Single Unit ',num2str(i)])
                    if i==size(sub(s).block(blocknum).bpseudo(lead).su,2)
                        xlabel('Time (ms)')
                        %if j==length(events)
                        %    leg=legend(err([1,4,5]),'Go','Turn','Swerve');
                        %    set(leg,'box','off','Position',[0.47,-0.01,0.1,0.05],'units','normalized','orientation','horizontal')
                        %    if size(sub(s).block(blocknum).bpseudo(lead).su,2)>1
                        %        %tightfig;
                        %    end
                        %end
                    else
                        %set(gca,'XTickLabel',{})
                    end
                    if j==1
                        ylabel('FR (Hz)')
                    else
                        set(gca,'YTickLabel',{})
                    end
                    set(gca,'box','off','FontSize',12)
                    %if i~=1
                    %
                    %end
                end
                ylims(1)=0.8*min(miny);
                ylims(2)=1.2*max(maxy);
                for j=1:length(events)
                    set(frplot(j),'YLim',ylims);
                end
            end
            suptitle(['Subject ',num2str(s),' Block ',num2str(blocknum),' Lead ',num2str(lead)])
            set(gcf,'Units','inches','Position',[3,3,12,3+2*size(sub(s).block(blocknum).bpseudo(lead).su,2)])
            print(['/Users/Dennis/Documents/MATLAB/Mogilner_Lab/interrupted_reach/Subject_',num2str(s),'_Block_',num2str(blocknum),'_Lead_',num2str(lead)],'-depsc')
        end
    end
end

%% generate matrix that will be input into PCA

eventsdpca=[1,2,3,9,4,8,7];
eventmap=[1,2,3,4,5,6,7];
types=[1,4,5];

eventsdpca_error=[1,2,3,9,7];
eventmap_error=[1,2,3,4,5];
types_error=[1,2];

%eventmap=[1,2,3,4,6,7];
%eventsdpca=[1,2,3,9,8,7];
%types=[4,5];
subevents=ismember(events,eventsdpca);
subevents_error=ismember(events,eventsdpca_error);
event_names_dpca=event_names(subevents);

for k=1:length(types)
    pretimessort(:,:,k)=cellfun(@(x,y) sort(x(y==types(k))),preall,typeallall,'uniformoutput',0);
    ntimes(:,:,k)=cellfun(@(x,y) sum(~isnan(x(y==types(k),:))),preall,typeallall,'uniformoutput',0);
end
pretimes=NaN*ones(size(pretimessort));
for k=1:length(types)
    for x=1:size(ntimes,1)
        for y=1:size(ntimes,2)
            if ntimes{x,y,k}~=0
                pretimes(x,y,k)=pretimessort{x,y,k}(round(ntimes{x,y,k}/3));
            end
        end
    end
end
pretimes=max(max(pretimes,[],3));

for k=1:length(types)
    posttimessort(:,:,k)=cellfun(@(x,y) sort(x(y==types(k))),postall,typeallall,'uniformoutput',0);
    ntimes(:,:,k)=cellfun(@(x,y) sum(~isnan(x(y==types(k),:))),postall,typeallall,'uniformoutput',0);
end
posttimes=NaN*ones(size(posttimessort));
for k=1:length(types)
    for x=1:size(ntimes,1)
        for y=1:size(ntimes,2)
            if ntimes{x,y,k}~=0
                posttimes(x,y,k)=posttimessort{x,y,k}(round(2*ntimes{x,y,k}/3));
            end
        end
    end
end
posttimes=min(min(posttimes,[],3));

pretimes(1:3)=[-0.5,-0.2,-0.2];
posttimes([1,2,7])=[0.3,0.3,1];

pretimes=pretimes(subevents);
posttimes=posttimes(subevents);

pretimes_error=pretimes(subevents_error);
posttimes_error=posttimes(subevents_error);

alltimes=sub(1).block(1).b(1).su(1,1).fr_times{1};
preinds=discretize(pretimes,alltimes);
postinds=discretize(posttimes,alltimes);

cuminds=[0,cumsum(postinds-preinds+1)];

preinds_error=discretize(pretimes_error,alltimes);
postinds_error=discretize(posttimes_error,alltimes);

cuminds_error=[0,cumsum(postinds_error-preinds_error+1)];


frs=zeros(1,2,sum(postinds-preinds+1));
frs_error=zeros(1,2,sum(postinds_error-preinds_error+1));
trialnum=zeros(1,2);
trialnum_error=trialnum;
count=0;
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            for i=1:size(sub(s).block(blocknum).bpseudo(lead).su,2)
                for k=1:length(types)
                    for j=1:length(eventsdpca)
                        
                        frs(count+1,k,cuminds(j)+1:cuminds(j+1))=sub(s).block(blocknum).bpseudo(lead).su(eventsdpca(j),i).fr{types(k)}(preinds(j):postinds(j));
                        trialnum(count+1,k)=size(sub(s).block(blocknum).bpseudo(lead).su(eventsdpca(j),i).fr_ind{types(k)},1);
                    end
                    
                    if k==3 || s==4
                        continue
                    end
                    for j=1:length(eventsdpca_error)
                         
                        frs_error(count+1,k,cuminds_error(j)+1:cuminds_error(j+1))=sub(s).block(blocknum).bpseudo(lead).su_error(eventsdpca_error(j),i).fr{types_error(k)}(preinds_error(j):postinds_error(j));
                        trialnum_error(count+1,k)=size(sub(s).block(blocknum).bpseudo(lead).su_error(eventsdpca_error(j),i).fr_ind{types_error(k)},1);
                    end
                    
                end
                s
                blocknum
                lead
                i
                count
                count=count+1;
            end
        end
    end
end

frtrial=NaN*ones(size(frs,1),size(frs,2),size(frs,3),max(max(trialnum)));
frtrial_error=NaN*ones(size(frs_error,1),size(frs_error,2),size(frs_error,3),max(max(trialnum_error)));

count=0;
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            for i=1:size(sub(s).block(blocknum).bpseudo(lead).su,2)
                for k=1:length(types)
                    origind=sub(s).block(blocknum).bpseudo(lead).su(1,i).trials(~fixbreakall{s,blocknum,lead,1,i}' & typeall{s,blocknum,lead,1,i}==types(k));
                    for j=1:length(eventsdpca)
                        inds=false(1,size(frtrial,4));
                        inds(1:length(origind))=ismember(origind,sub(s).block(blocknum).bpseudo(lead).su(eventsdpca(j),i).trials(~fixbreakall{s,blocknum,lead,eventmap(j),i}' & typeall{s,blocknum,lead,eventmap(j),i}==types(k)));
                        frtrial(count+1,k,cuminds(j)+1:cuminds(j+1),inds)=permute(sub(s).block(blocknum).bpseudo(lead).su(eventsdpca(j),i).fr_ind{types(k)}(:,preinds(j):postinds(j))',[3,4,1,2]);
                    end
                    
                    if k==3 || s==4
                        continue
                    end
                    origind=sub(s).block(blocknum).bpseudo(lead).su(1,i).trials(~fixbreakall_error{s,blocknum,lead,1,i}' & typeall_error{s,blocknum,lead,1,i}==types_error(k));
                    for j=1:length(eventsdpca_error)
                        inds=false(1,size(frtrial_error,4));
                        inds(1:length(origind))=ismember(origind,sub(s).block(blocknum).bpseudo(lead).su(eventsdpca_error(j),i).trials(~fixbreakall_error{s,blocknum,lead,eventmap_error(j),i}' & typeall_error{s,blocknum,lead,eventmap_error(j),i}==types_error(k)));
                        frtrial_error(count+1,k,cuminds_error(j)+1:cuminds_error(j+1),inds)=permute(sub(s).block(blocknum).bpseudo(lead).su_error(eventsdpca_error(j),i).fr_ind{types_error(k)}(:,preinds_error(j):postinds_error(j))',[3,4,1,2]);
                    end
                end
                count=count+1;
            end
        end
    end
end

%exclude sessions with too few trials
frs=frs([1:11,13:24,27,28,33:end],:,:);
frtrial=frtrial([1:11,13:24,27,28,33:end],:,:,:);
trialnum=trialnum([1:11,13:24,27,28,33:end],:);

frs_error=frs_error([1:11,13:15,17:24,27,28,33:end],:,:);
frtrial_error=frtrial_error([1:11,13:15,17:24,27,28,33:end],:,:,:);
trialnum_error=trialnum_error([1:11,13:15,17:24,27,28,33:end],:);

%frtrial=bsxfun(@times,frtrial,1./mean(mean(frs,3),2));
%frs=bsxfun(@times,frs,1./mean(mean(frs,3),2));

combinedParams={{1,[1,2]},{2}};
%notToSplit={{3,[2,3]},{[1,3],[1,2,3]}};
margNames={'Trial Types','Condition-Independent'};

time=0:0.01:0.01*(length(frs)-1);
%timeEvents=time(cuminds(2:end-1));
timeEvents=time((201-preinds)+cuminds(1:end-1));

frmerged=reshape(permute(frtrial,[1,3,4,2]),35,353,[],1);
mu=mean(reshape(frmerged(:,:,:),35,[]),2,'omitnan');
stdev=std(reshape(frmerged(:,:,:),35,[]),[],2,'omitnan');
frtrial=bsxfun(@times,bsxfun(@minus,frtrial,mu),1./stdev);
frs=squeeze(mean(frtrial,4,'omitnan'));

%%
dataDim=size(frs);
Xcen=frs(:,:);
Xcen=bsxfun(@minus, Xcen, mean(Xcen,2));
[coeff,score,~,~,explained]=pca(Xcen');

proj1=reshape(score(:,1),dataDim(2:end));
proj2=reshape(score(:,2),dataDim(2:end));
proj3=reshape(score(:,3),dataDim(2:end));

for ii=1:3
    subplot(3,2,(ii-1)*2+1)
   
    
    for j=1:7
        if mod(j,2)==1
            patch([time(cuminds(j)+1),time(cuminds(j+1)),time(cuminds(j+1)),time(cuminds(j)+1)],...
                [-2,-2,6,6],[0.8,0.8,0.8],'FaceAlpha',0.5,'EdgeColor','none')
            hold on
        end
    end
    
    l((ii-1)*3+(1:3))=plot(time,reshape(score(:,ii),dataDim(2:end)),'LineWidth',2);
    ylim([-2,6])
    xlim([min(time),max(time)])
    set(gca,'box','off','FontSize',12)
    ylabel({'Z-Scored'; 'Firing Rate'})
    if ii==2
        xlabel('Time (s)')
    elseif ii==1
        text(0.5,1.3,'Top 3 Principal Components','FontSize',16,'units','normalized','horizontalAlignment','center')
    end
    %if ii==2
    %    leg1=legend(l(4:6),'Reach','Planned Turn','Impromptu Turn');
    %    set(leg1,'box','off')
    %    leg1.Position(1:2)=[0.4,0.55];
    
    %end
    
    line(repmat(timeEvents,2,1),repmat([-2;3.5],1,length(timeEvents)),'color',[0,0,0.5])
    for ev=1:length(timeEvents)
        if ev<5
            text(timeEvents(ev)-0.05,4,event_names{eventmap(ev)},'horizontalAlignment','center','rotation',85,'Color',[0,0,0.5],'FontSize',8)
        else
            text(timeEvents(ev)+0.05,4,event_names{eventmap(ev)},'horizontalAlignment','center','rotation',275,'Color',[0,0,0.5],'FontSize',8)
        end
    end
    
    if ii==1
        
        plot(time(h),-1*ones(nnz(h)),'k','LineWidth',3)
        
    end
        
        
end
%leg1.String=leg1.String(1:3);

subplot(3,2,2)
plot(cumsum(explained),'.-','LineWidth',2,'MarkerSize',10,'color','k')
ylim([0 100])
xlabel('Number of Principal Components')
ylabel({'% Cumulative'; 'Explained Variance'})
set(gca,'box','off','FontSize',12)


ax=subplot(3,2,4);
cdata=arrayfun(@(x,y) linspace(x,y,353),mean(cat(3,colors,ones(size(colors))),3),...
    mean(cat(3,colors,zeros(size(colors))),3),'uniformoutput',0);
cdata=reshape(vertcat(cdata{:}),7,3,353);

markers={'v','o','^','h','s','o','p'};
markersize([3,5,6,7])=[13,13,13,15];
evinds=cuminds(2:end)-(postinds-200);
for ii=1:3
surface([proj1(ii,:);proj1(ii,:)],[proj2(ii,:);proj2(ii,:)],zeros(2,353),...
    permute(cat(1,cdata(ii,:,:),cdata(ii,:,:)),[1,3,2]),'facecol','no','edgecol','interp','linew',2,'Marker','.','MarkerSize',10);
end
hold on
for ii=1:3
    for j=[3,5,6,7]
    plot(proj1(ii,evinds(j)),proj2(ii,evinds(j)),'marker',markers{j},'markeredgecolor','k','markerfacecolor',colors(ii,:),'markersize',markersize(j),'linestyle','none','LineWidth',2)
    end
end
text(0.5,-0.1,'PC 1','units','normalized','horizontalalignment','center','FontSize',12)
text(-0.1,0.5,'PC 2','units','normalized','horizontalalignment','center','FontSize',12,'rotation',90)
xlim([min(min(proj1)),max(max(proj1))])
ylim([min(min(proj2)),max(max(proj2))])
colormap(ax,'gray')
cbar=colorbar('north','Ticks',[],'box','off','direction','reverse');
cbar.Position=[0.8,0.6,0.1,0.01];
annotation('textarrow',[0.85,0.9],[0.62,0.62],'String','Time  ','headstyle','cback3')
set(gca,'box','off','FontSize',12,'Visible','off')


subplot(3,2,6)
for ii=1:3
surface([proj1(ii,:);proj1(ii,:)],[proj3(ii,:);proj3(ii,:)],zeros(2,353),...
    permute(cat(1,cdata(ii,:,:),cdata(ii,:,:)),[1,3,2]),'facecol','no','edgecol','interp','linew',2,'Marker','.','MarkerSize',10);
end
hold on
for ii=1:3
    for j=[3,5,6,7]
    plot(proj1(ii,evinds(j)),proj3(ii,evinds(j)),'marker',markers{j},'markeredgecolor','k','markerfacecolor',colors(ii,:),'markersize',markersize(j),'linestyle','none','LineWidth',2)
    m(j)=plot(0,0,'marker',markers{j},'markeredgecolor','k','markerfacecolor','none','markersize',markersize(j),'linestyle','none','LineWidth',2,'visible','off');
    end
end
text(0.5,-0.1,'PC 1','units','normalized','horizontalalignment','center','FontSize',12)
text(-0.1,0.5,'PC 3','units','normalized','horizontalalignment','center','FontSize',12,'rotation',90)
xlim([min(min(proj1)),max(max(proj1))])
ylim([min(min(proj3)),max(max(proj3))])
leg=legend(m([3,5,6,7]),'\color{black}Go','\color{black}Turn Cue','\color{black}Turn Onset','\color{black}Feedback');
set(leg,'box','off','units','normalized');
leg.Position(1:2)=[0.6,0.3];
set(gca,'box','off','FontSize',12,'Visible','off')

%% plot local LFPs aligned to events
events=[1:3,8];
type=[1,4,5];
t=-2:0.001:2;
event_names={'Fixation','Target','Go','Turn'};
figure
for s=[1,2,4,7]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            figure
            for j=1:length(events)
                trialtype=sub(s).block(blocknum).b(lead).type;
                trialtype=trialtype(sub(s).block(blocknum).b(lead).su(events(j),1).trials);
                subplot(length(events),1,j)
                if ismember(j,1:3)
                    col=1;
                    for k=1:3
                        meansig=mean(cell2mat(sub(s).block(blocknum).ao(lead).trials(events(j),trialtype==type(k))'));
                        se=std(cell2mat(sub(s).block(blocknum).ao(lead).trials(events(j),trialtype==type(k))'))./sqrt(nnz(trialtype==type(k)));
                        
                        meansig=meansig(8000:12000);
                        se=se(8000:12000);
                        
                        plot(t,meansig,'color',colors(k,:))
                        hold on
                        h(k)=patch([t,fliplr(t)],[meansig-se,fliplr(meansig+se)],colors(col,:),'FaceAlpha',0.5,'EdgeColor','none');
                        col=col+1;
                    end
                elseif j==4
                    col=2;
                    for k=2:3
                        meansig=mean(cell2mat(sub(s).block(blocknum).ao(lead).trials(events(j),trialtype==type(k))'));
                        se=std(cell2mat(sub(s).block(blocknum).ao(lead).trials(events(j),trialtype==type(k))'))./sqrt(nnz(trialtype==type(k)));
                        
                        meansig=meansig(8000:12000);
                        se=se(8000:12000);
                        
                        plot(t,meansig,'color',colors(k,:))
                        hold on
                        patch([t,fliplr(t)],[meansig-se,fliplr(meansig+se)],colors(col,:),'FaceAlpha',0.5,'EdgeColor','none')
                        col=col+1;
                    end
                    xlabel('Time (s)')
                end
                ylim([-0.3,0.1])
                xlim([-2,2])
                title(event_names(j))
                set(gca,'box','off')
                %plot(meansig)
            end
        end
    end
end
leg=legend(h,'Go','Turn','Swerve');
set(leg,'Position',[0.35,-0.01,0.3,0.1],'box','off','orientation','horizontal')

%% LFP

%% calculate continuous PSTHs separated by trial type (requires chronux toolbox)
event_ord=[1,2,3,10,4,9,8];
suevent_ord=[1,2,3,9,4,8,7];
types=[1,4,5];

preadj=[0.5,0.1,0.1,0.1,0.1,0.1,0.1];
postadj=[0,0,-0.15,0,-0.15,0,0];

preall=cell(11,1);
postall=preall;

presplit=cell(11,3);
postsplit=presplit;
typeall=cell(11,1);
fixbreakall=presplit;

for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        result=sub(s).block(blocknum).b(1).result;
        type=sub(s).block(blocknum).b(1).type;
        trialtimes=sub(s).block(blocknum).b(1).t;
        tev1_trials=sub(s).block(blocknum).b(1).tev1_trials;
        tev1_trials(type==4,4)=NaN;
        event_happened=sub(s).block(blocknum).b(1).event_happened;
        event_happened=[event_happened,tev1_trials(:,9)~=-1 & ~isnan(tev1_trials(:,9)),tev1_trials(:,10)~=-1 & ~isnan(tev1_trials(:,10))];
        event_happened(type==4,4)=0;
        event_happened=event_happened(:,event_ord);
        tev1_trials=tev1_trials(:,event_ord);
        result(result==2 & type==4 & ~(tev1_trials(:,6)==-1 | isnan(tev1_trials(:,6)))')=1; %turn trials where a turn was identified but were marked incorrect
        result(result==2 & type==5 & ~isnan(tev1_trials(:,5))' & ~(tev1_trials(:,6)==-1 | isnan(tev1_trials(:,6)))')=1; %swerve trials where a swerve signal was present and a turn was identified but were marked incorrect
        result(result==1 & (type==4 | type==5) & tev1_trials(:,6)'==-1)=2; %turn or swerve trials marked correct but no turn identified
        result(~isnan(tev1_trials(:,6)) & tev1_trials(:,6)~=-1 & tev1_trials(:,6)<tev1_trials(:,3))=4; %trials where a turn was identified before the go cue (i.e. incorrectly identified), ideally would go back and fix turn identification here
        tev1_trials(result'==1 & tev1_trials(:,4)<tev1_trials(:,3),4)=tev1_trials(result'==1 & tev1_trials(:,4)<tev1_trials(:,3),3)+0.05;
        %fixbreak=tev1_trials(:,4)<tev1_trials(:,3) | isnan(tev1_trials(:,2)) | isnan(tev1_trials(:,3));
        %error=tev1_trials(:,6)==-1 | (type==1 & result~=1)';
        %fixbreak=fixbreak|error;
        fixbreak=result'==2|result'==3|result'==4;
        %fixbreak=false(size(fixbreak));
        
        pre=NaN*ones(size(tev1_trials));
        post=pre;
        
        
        for ev=1:length(event_ord)
            j=event_ord(ev);
            switch j
                case 1
                    pre(1,ev)=-5;
                    pre(2:end,ev)=trialtimes(1:end-1,2)-trialtimes(2:end,1)+preadj(ev);
                    
                    [row,col]=find(event_happened(:,ev+1:end));
                    nextev=accumarray(row,col,[size(event_happened,1),1],@min,NaN);
                    nextev=nextev+ev;
                    post(:,ev)=tev1_trials(sub2ind(size(tev1_trials),1:length(type),nextev'))'+postadj(ev)-tev1_trials(:,ev);
                    
                case 8
                    [row,col]=find(event_happened(:,1:ev-1));
                    priorev=accumarray(row,col,[size(event_happened,1),1],@max,NaN);
                    pre(:,ev)=tev1_trials(sub2ind(size(tev1_trials),1:length(type),priorev'))'+preadj(ev)-tev1_trials(:,ev);
                    
                    post(1:end-1,ev)=trialtimes(2:end,1)-trialtimes(1:end-1,2)+postadj(ev);
                    post(end,ev)=5;
                otherwise
                    [row,col]=find(event_happened(:,1:ev-1));
                    priorev=accumarray(row,col,[size(event_happened,1),1],@max,NaN);
                    pre(:,ev)=tev1_trials(sub2ind(size(tev1_trials),1:length(type),priorev'))'+preadj(ev)-tev1_trials(:,ev);
                    
                    [row,col]=find(event_happened(:,ev+1:end));
                    nextev=accumarray(row,col,[size(event_happened,1),1],@min,NaN);
                    nextev=nextev+ev;
                    post(:,ev)=tev1_trials(sub2ind(size(tev1_trials),1:length(type),nextev'))'+postadj(ev)-tev1_trials(:,ev);
                    
            end
            presplit{s,blocknum}=pre;
            postsplit{s,blocknum}=post;
            %fixbreakall{s,blocknum}=fixbreak;
            %typeall{s,blocknum}=type;
            
            pre(logical(~event_happened))=NaN;
            post(logical(~event_happened))=NaN;
            
            preall{s}=[preall{s};pre];
            postall{s}=[postall{s};post];
            
            fixbreakall{s}=[fixbreakall{s};fixbreak];
            typeall{s}=[typeall{s};type'];
            
        end
        
        sub(s).block(blocknum).presplit=pre;
        sub(s).block(blocknum).postsplit=post;
        sub(s).block(blocknum).fixbreak=fixbreak;
        sub(s).block(blocknum).type=type;
        sub(s).block(blocknum).tev1_trials=tev1_trials;
        
        
    end
end
%%
tic
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            mfile=matfile([pathbase,'/LFP/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(lead)]);
            lfp=psthc(sub(1).pathbase,sub(s).block(blocknum),mfile,2000,[-2,2],event_ord,'LFP_align','lfp',logical(sub(s).block(blocknum).ao.artlog(lead,:)));
            toc
        end
    end
end
%%
tic
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            mfile=matfile([pathbase,'/Hilbert/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(lead)]);
            hil=psthc(sub(1).pathbase,sub(s).block(blocknum),mfile,2000,[-2,2],event_ord,'Hilbert_align','hil',logical(sub(s).block(blocknum).ao.artlog(lead,:)));
            toc
        end
    end
end
%%
tic
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            s
            blocknum
            lead
            mfile=matfile([pathbase,'/CWT/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(lead)]);
            wt=psthc(sub(1).pathbase,sub(s).block(blocknum),mfile,2000,[-2,2],event_ord,'CWT_align','wt',logical(sub(s).block(blocknum).ao.artlog(lead,:)));
            toc
        end
    end
end
%% plot lfp
fs=2000;
T=[-2,2];
colorsall=get(gca,'colororder');
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        %figure
        pre=sub(s).block(blocknum).presplit;
        post=sub(s).block(blocknum).postsplit;
        fixbreak=sub(s).block(blocknum).fixbreak;
        type=sub(s).block(blocknum).type;
        
        bounds=[median(pre(~fixbreak,:),'omitnan');median(post(~fixbreak,:),'omitnan')];
        bounds(bounds<T(1))=T(1);
        bounds(bounds>T(2))=T(2);
        bound_disc=fix(bounds*fs)+fs/2*(T(2)-T(1))+1;
        bounds=fix(bounds*fs)/fs;
        for  lead=1:length(sub(s).block(blocknum).b)
            mfile=matfile([pathbase,'/LFP_align/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(lead)]);
            yrange=[0,0];
            lfpmean=NaN*ones(3,7,8001);
            lfperrmat=lfpmean;
            lfpbase=squeeze(mfile.lfp(:,1,preinds(1):preinds(1)+1000));
            lfpmu=mean(lfpbase(:),'omitnan');
            lfpstd=std(lfpbase(:),'omitnan');
            for j=1:length(events)
                %h(j)=subplot(length(sub(s).block(blocknum).b),length(events),(lead-1)*length(events)+j);
                lfpfull=squeeze(mfile.lfp(:,j,bound_disc(1,j):bound_disc(2,j)));
                lfpfull=(lfpfull-lfpmu)./lfpstd;
                event=events(j);
                switch event
                    case 4
                        types=5;
                        colors=colorsall(3,:);
                    case 8
                        types=[4,5];
                        colors=colorsall(2:3,:);
                    otherwise
                        types=[1,4,5];
                        colors=colorsall(1:3,:);
                        
                end
                for k=1:length(types)
                    lfp=squeeze(mean(lfpfull(type(~fixbreak)==types(k),:),'omitnan'));
                    lfpmean(k,j,bound_disc(1,j):bound_disc(2,j))=lfp;
                    lfperr=squeeze(std(lfpfull(type(~fixbreak)==types(k),:),'omitnan')./sqrt(sum(~isnan(lfpfull(type(~fixbreak)==types(k),:)))));
                    lfperrmat(k,j,bound_disc(1,j):bound_disc(2,j))=lfperr;
                    %patch([bounds(1,j):1/fs:bounds(2,j),fliplr(bounds(1,j):1/fs:bounds(2,j))],...
                    %    [lfp-lfperr,fliplr(lfp+lfperr)],mean([1,1,1;colors(k,:)]),'EdgeColor','none','FaceAlpha',0.5);
                    %hold on
                    %plot(bounds(1,j):1/fs:bounds(2,j),lfp,'color',colors(k,:))
                    %yrange(1)=min([yrange(1),lfp]);
                    %yrange(2)=max([yrange(2),lfp]);
                end
                
            end
            %for j=1:length(events)
            %    set(h(j),'ylim',yrange);
            %end
            save([pathbase,'/LFP_mean/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(lead)],'lfpmean','lfperrmat','-v7.3');
        end
        %suptitle(['Subject ',num2str(s),' Block ',num2str(blocknum)]);
        
    end
end

%%
fs=2000;
T=[-2,2];
colorsall=get(gca,'colororder');
tic
for s=[1,2,4,7,8,9,11]
    s
    for blocknum=1:length(sub(s).block)
        blocknum
        pre=sub(s).block(blocknum).presplit;
        post=sub(s).block(blocknum).postsplit;
        fixbreak=sub(s).block(blocknum).fixbreak;
        type=sub(s).block(blocknum).type;
        ntypes=zeros(1,3);
        for a=1:3
            ntypes(a)=nnz(type(~fixbreak)==types(a));
        end
        maxn=max(ntypes);
        summaxn=sum(maxn);
        nreps=round(exp(-(1:maxn)/20)*50+25);
        
        
        bounds=[median(pre(~fixbreak,:),'omitnan');median(post(~fixbreak,:),'omitnan')];
        bounds(bounds<T(1))=T(1);
        bounds(bounds>T(2))=T(2);
        bound_disc=fix(bounds*fs)+fs/2*(T(2)-T(1))+1;
        bounds=fix(bounds*fs)/fs;
        for  lead=1:length(sub(s).block(blocknum).b)
            lead
            %figure
            mfile=matfile([pathbase,'/CWT_align/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(lead)]);
            yrange=[0,0];
            wtfull=mfile.wt;
            
            wtbaseevoked=NaN*ones(55,maxn,max(nreps));
            wtbaseinduced=wtbaseevoked;
            for ii=1:maxn
                for jj=1:nreps(ii)
                    wtbaseevoked(:,ii,jj)=mean(abs(mean(wtfull(randperm(summaxn,ii),1,bound_disc(1,1):(bound_disc(1,1)+500),:),'omitnan')));
                    wtbaseinduced(:,ii,jj)=mean(mean(abs(wtfull(randperm(summaxn,ii),1,bound_disc(1,1):(bound_disc(1,1)+500),:)),'omitnan'));
                end
            end
            wtbasepower=wtbaseinduced;
            wtbaseinduced=mean(wtbaseinduced-wtbaseevoked,3,'omitnan');
            wtbaseevoked=mean(wtbaseevoked,3,'omitnan');
            
            %wtbaseevoked=squeeze(mean(abs(mean(wtfull(:,1,bound_disc(1,1):(bound_disc(1,1)+500),:),'omitnan'))));
            %wtbaseinduced=squeeze(mean(mean(abs(wtfull(:,1,bound_disc(1,1):(bound_disc(1,1)+500),:)),'omitnan')));
            wtevoked=NaN*ones(3,7,8001,55);
            wtinduced=wtevoked;
            wtpower=wtevoked;
            for j=1:length(events)
                event=events(j);
                %wtfull=squeeze(mfile.wt(:,j,bound_disc(1,j):bound_disc(2,j),:));
                switch event
                    case 4
                        types=5;
                        colors=colorsall(3,:);
                        plotrow=3;
                    case 8
                        types=[4,5];
                        colors=colorsall(2:3,:);
                        plotrow=[2,3];
                    otherwise
                        types=[1,4,5];
                        colors=colorsall(1:3,:);
                        plotrow=1:3;
                        
                end
                for k=1:length(types)
                    %h(j,k)=subplot(3,length(events),(plotrow(k)-1)*length(events)+j);
                    bounds=[median(pre(~fixbreak & type'==types(k),:),'omitnan');median(post(~fixbreak & type'==types(k),:),'omitnan')];
                    bounds(bounds<T(1))=T(1);
                    bounds(bounds>T(2))=T(2);
                    bound_disc=fix(bounds*fs)+fs/2*(T(2)-T(1))+1;
                    bounds=fix(bounds*fs)/fs;
                    evok=squeeze(mean(wtfull(type(~fixbreak)==types(k),j,bound_disc(1,j):bound_disc(2,j),:),'omitnan'));
                    evok=abs(evok)';
                    evoknnan=squeeze(sum(~isnan(wtfull(type(~fixbreak)==types(k),j,bound_disc(1,j):bound_disc(2,j),1))));
                    [~,maxind]=max(sum(squeeze(sum(~isnan(wtfull(type(~fixbreak)==types(k),j,bound_disc(1,j):bound_disc(2,j),:)))),2));
                    correction=median(bsxfun(@times,squeeze(evok(20:55,:)),1./evok(20:55,maxind)));
                    %evok=bsxfun(@times,evok,1./correction);
                    induc=squeeze(mean(abs(wtfull(type(~fixbreak)==types(k),j,bound_disc(1,j):bound_disc(2,j),:)),'omitnan'))';
                    %induc=(bsxfun(@times,induc'-evok,1./wtbaseinduced));
                    %induc(:,evoknnan>0)=(induc(:,evoknnan>0)-evok(:,evoknnan>0))./bsxfun(@minus,wtbaseinduced,wtbaseevoked(:,evoknnan(evoknnan>0)));
                    pow=induc;
                    pow(:,evoknnan>0)=pow(:,evoknnan>0)/wtbasepower(:,evoknnan(evoknnan>0));
                    induc(:,evoknnan>0)=(induc(:,evoknnan>0)-evok(:,evoknnan>0))./wtbaseinduced(:,evoknnan(evoknnan>0));
                    %evok=(bsxfun(@times,evok,1./wtbaseevoked));
                    evok(:,evoknnan>0)=evok(:,evoknnan>0)./wtbaseevoked(:,evoknnan(evoknnan>0));
                    wtevoked(k,j,bound_disc(1,j):bound_disc(2,j),:)=evok';
                    wtinduced(k,j,bound_disc(1,j):bound_disc(2,j),:)=induc';
                    wtpower(k,j,bound_disc(1,j):bound_disc(2,j),:)=pow';
                    %surf(bounds(1,j):1/fs:bounds(2,j),f,log10(evok),'edgecolor','none')
                    %view(0,90);
                    %set(gca,'yscale','log','YTick',[2:2:10,20:20:60,100,200,400])
                    %xlim(bounds(:,j)')
                    %hold on
                    %yrange(1)=min([yrange(1);log10(evok(:))]);
                    %yrange(2)=max([yrange(2);log10(evok(:))]);
                end
                
            end
            toc
            %yrange=[-1,1];
            %for j=1:length(events)
            %    event=events(j);
            %    switch event
            %        case 4
            %            plotrow=3;
            %        case 8
            %            plotrow=[2,3];
            %        otherwise
            %            plotrow=1:3;
            %    end
            %    for k=1:length(plotrow)
            %        caxis(h(j,k),yrange);
            %    end
            %end
            %suptitle(['Subject ',num2str(s),' Block ',num2str(blocknum),' Lead ',num2str(lead)]);
            save([pathbase,'/CWT_mean/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(lead)],'wtevoked','wtinduced','-v7.3');
        end
    end
end

%%
fs=2000;
T=[-2,2];
colorsall=get(gca,'colororder');

for k=1:3
    pretimessort(:,k)=cellfun(@(x,y) sort(x(y==types(k),:)),preall,typeall,'uniformoutput',0);
    ntimes(:,k)=cellfun(@(x,y) sum(~isnan(x(y==types(k),:))),preall,typeall,'uniformoutput',0);
end
pretimes=NaN*ones(length(pretimessort),7,3);
for k=1:3
    for x=1:length(ntimes)
        for y=1:length(ntimes{x,k})
            if ntimes{x,k}(y)~=0
                pretimes(x,y,k)=pretimessort{x,k}(round(ntimes{x,k}(y)/3),y);
            end
        end
    end
end

pretimes=max(max(pretimes,[],3));

for k=1:3
    posttimessort(:,k)=cellfun(@(x,y) sort(x(y==types(k),:)),postall,typeall,'uniformoutput',0);
    ntimes(:,k)=cellfun(@(x,y) sum(~isnan(x(y==types(k),:))),postall,typeall,'uniformoutput',0);
end
posttimes=NaN*ones(length(posttimessort),7,3);
for k=1:3
    for x=1:length(ntimes)
        for y=1:length(ntimes{x,k})
            if ntimes{x,k}(y)~=0
                posttimes(x,y,k)=posttimessort{x,k}(round(2*ntimes{x,k}(y)/3),y);
            end
        end
    end
end

posttimes=min(min(posttimes,[],3));

alltimes=T(1):1/fs:T(2);
preinds=discretize(pretimes,alltimes);
postinds=discretize(posttimes,alltimes);

wtevoked=NaN*ones(3,7,8001,55,22);
wtinduced=wtevoked;
tic
count=1;
for s=[1,2,4,7,8,9,11]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            [s,blocknum,lead,count]
            mfile=matfile([pathbase,'/CWT_mean/Sub',num2str(s),'_Block',num2str(blocknum),'_Lead',num2str(lead)]);
            for j=1:7
                wtevoked(:,j,preinds(j):postinds(j),:,count)=mfile.wtevoked(:,j,preinds(j):postinds(j),:);
                wtinduced(:,j,preinds(j):postinds(j),:,count)=mfile.wtinduced(:,j,preinds(j):postinds(j),:);
            end
            count=count+1;
            toc
        end
    end
end

wtevokedmean=squeeze((mean(wtevoked(:,:,:,:,[1:5,7:13,15:18:end]),5,'omitnan')));
wtinducedmean=squeeze((mean(wtinduced(:,:,:,:,[1:5,7:13,15:18:end]),5,'omitnan')));

%%
trial_name={'Go Trials','Turn Trials','Swerve Trials'};
yrange=[0,0];
tic
%for ii=1:22
figure
for j=1:length(events)
    event=events(j);
    switch event
        case 4
            types=5;
            colors=colorsall(3,:);
            plotrow=3;
        case 8
            types=[4,5];
            colors=colorsall(2:3,:);
            plotrow=[2,3];
        otherwise
            types=[1,4,5];
            colors=colorsall(1:3,:);
            plotrow=1:3;
            
    end
    for k=1:length(types)
        h(j,k)=subplot(3,length(events),(plotrow(k)-1)*length(events)+j);
        surf(alltimes(preinds(j):postinds(j)),f,squeeze((wtevokedmean(k,j,preinds(j):postinds(j),:)))','edgecolor','none')
        hold on
        l=line([0,0],[2,500],[2,2]);
        set(l,'LineWidth',2,'color','k');
        view(0,90);
        set(gca,'yscale','log','YTick',[2,4,10,30,100,300],'FontSize',14)
        ylim([2,500])
        if types(k)~=5
            %set(gca,'XTickLabel',{});
        else
            xlabel('Time (s)')
            text(0.5,3.85,event_names{j},'units','normalized','horizontalAlignment','center','FontWeight','bold','FontSize',14)
        end
        
        if j~=1
            set(gca,'YTickLabel',{});
        else
            ylabel('Frequency (Hz)','FontSize',12)
            text(-1,0.5,trial_name{k},'rotation',90,'FontSize',18,...
                'FontWeight','bold','units','normalized','horizontalAlignment','center',...
                'Color',colors(k,:))
        end
        if j==7 && k==1
            cbar=colorbar('WestOutside');
            cbar.Position(1:2)=[0.63,0.4];
            cbarlabel=annotation('textarrow',[0.65,0.65],[0.397,0.397],'String',...
                'Normalized Induced Power','units','normalized','LineStyle','none',...
                'HeadStyle','none','TextRotation',270,'FontSize',10);
        end
        
        xlim([pretimes(j),posttimes(j)])
        hold on
        %yrange(1)=min([yrange(1);(wtevoked(:))]);
        %yrange(2)=max([yrange(2);(wtevoked(:))]);
        toc
    end
    
    
end

for j=1:length(events)
    event=events(j);
    switch event
        case 4
            plotrow=3;
        case 8
            plotrow=[2,3];
        otherwise
            plotrow=1:3;
    end
    for k=1:length(plotrow)
        %caxis(h(j,k),yrange);
        caxis(h(j,k),[0,2.25]);
    end
end
%end
