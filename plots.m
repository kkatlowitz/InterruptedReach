%% calculate PSTHs separated by trial type (requires chronux toolbox)
for s=[1,2,4,7]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            for j=1:8
                for i=1:size(sub(s).block(blocknum).b(lead).su,2)
                    
                    type=sub(s).block(blocknum).b(lead).type;
                    type=type(sub(s).block(blocknum).b(lead).su(j,i).trials);
                    data=sub(s).block(blocknum).b(lead).su(j,i).data;
                    for k=[1,4,5]
                        if nnz(type==k)~=0
                            [sub(s).block(blocknum).b(lead).su(j,i).fr{k},sub(s).block(blocknum).b(lead).su(j,i).fr_times{k},...
                                sub(s).block(blocknum).b(lead).su(j,i).fr_se{k}]=...
                                psth(data(type==k),0.05,'n',[-4,4],2,-4:0.01:4);
                        end
                    end
                end
            end
            
        end
    end
end
%% plot and save rasters with PSTH aligned to particular event

figure
colors=get(gca,'colororder');
events=[1,2,3,8];
event_names={'Fixation','Target','Go','Turn'};

for s=[1,2,4,7]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            if ~(s==1 && blocknum==1 && lead==1)
                figure
            end
            
            for j=1:length(events)
                for i=1:size(sub(s).block(blocknum).b(lead).su,2)
                    if length(sub(s).block(blocknum).b(lead).su(events(j),i).strials)<100 % || ~isfield(sub(s).block(blocknum).b(lead).su(events(j),i),'fr')
                        continue
                    end
                    subplot(2*length(events),size(sub(s).block(blocknum).b(lead).su,2),2*(j-1)*size(sub(s).block(blocknum).b(lead).su,2)+i)
                    plot(sub(s).block(blocknum).b(lead).su(events(j),i).stimes_align*1000,sub(s).block(blocknum).b(lead).su(events(j),i).strials,'.','MarkerSize',1,'MarkerEdgeColor','k') %plot rasters
                    xlim([-2000,2000])
                    ylim([0,max(sub(s).block(blocknum).b(lead).su(events(j),i).strials)+1])
                    if j==1
                        title(['Single Unit ',num2str(i)])
                    end
                    
                    if i==1
                        text(-0.1,-0.25,event_names{j},'units','normalized','rotation',90,'Fontweight','bold','horizontalalignment','center')
                    else
                        set(gca,'YTickLabel',{})
                    end
                    set(gca,'box','off','XTickLabel',{});
                    subplot(2*length(events),size(sub(s).block(blocknum).b(lead).su,2),2*(j-1)*size(sub(s).block(blocknum).b(lead).su,2)+size(sub(s).block(blocknum).b(lead).su,2)+i)
                    
                    
                    if j==4
                        
                        
                        col=2;
                        for k=[4,5]
                            patchx=[sub(s).block(blocknum).b(lead).su(events(j),i).fr_times{k}*1000,fliplr(sub(s).block(blocknum).b(lead).su(events(j),i).fr_times{k}*1000)];
                            patchy=[sub(s).block(blocknum).b(lead).su(events(j),i).fr{k}-sub(s).block(blocknum).b(lead).su(events(j),i).fr_se{k},fliplr(sub(s).block(blocknum).b(lead).su(events(j),i).fr{k}+sub(s).block(blocknum).b(lead).su(events(j),i).fr_se{k})];
                            err(k)=patch(patchx(~isnan(patchy)),patchy(~isnan(patchy)),mean([colors(col,:);1,1,1]),'EdgeColor','none','FaceAlpha',0.5);
                            
                            hold on
                            
                            l(k)=plot(sub(s).block(blocknum).b(lead).su(events(j),i).fr_times{k}*1000,sub(s).block(blocknum).b(lead).su(events(j),i).fr{k},'Color',colors(col,:));
                            col=col+1;
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
                        end
                    end
                    xlim([-2000,2000])
                    ylims=get(gca,'YLim');
                    ylims(1)=0;
                    set(gca,'YLim',ylims);
                    if j==4
                        xlabel('Time (ms)')
                        if i==size(sub(s).block(blocknum).b(lead).su,2)
                            leg=legend(err([1,4,5]),'Go','Turn','Swerve');
                            set(leg,'box','off','Position',[0.47,-0.01,0.1,0.05],'units','normalized','orientation','horizontal')
                        end
                    else
                        set(gca,'XTickLabel',{})
                    end
                    set(gca,'box','off')
                end
            end
            suptitle(['Subject ',num2str(s),' Block ',num2str(blocknum),' Lead ',num2str(lead)])
            set(gcf,'Units','inches','Position',[3,3,6*size(sub(s).block(blocknum).b(lead).su,2),10])
            print(['/Users/Dennis/Documents/MATLAB/Mogilner_Lab/interrupted_reach/subject_',num2str(s),'_Block_',num2str(blocknum),'_Lead_',num2str(lead)],'-depsc')
        end
    end
end
%% plot PSTHs for all events
%s: sub
%blocknum: block
%lead: lead
%i: single unit
%j: event
%k: trial type

event_plot=[1:3,8];
types=[1,4,5];
event_names={'Fixation','Target','Go','Turn'};
figure
colors=get(gca,'colororder');
for s=[1,2,4,7]
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            figure
            for j=1:length(event_plot)
                for i=1:size(sub(s).block(blocknum).b(lead).su,2)
                    subplot(length(event_plot),size(sub(s).block(blocknum).b(lead).su,2),size(sub(s).block(blocknum).b(lead).su,2)*(j-1)+i)
                    
                    for k=1:length(types)
                        if ~isempty(sub(s).block(blocknum).b(lead).su(event_plot(j),i).fr_times{types(k)})
                            patchx=[sub(s).block(blocknum).b(lead).su(event_plot(j),i).fr_times{types(k)},fliplr(sub(s).block(blocknum).b(lead).su(event_plot(j),i).fr_times{types(k)})];
                            patchy=[sub(s).block(blocknum).b(lead).su(event_plot(j),i).fr{types(k)}-sub(s).block(blocknum).b(lead).su(event_plot(j),i).fr_se{types(k)},fliplr(sub(s).block(blocknum).b(lead).su(event_plot(j),i).fr{types(k)}+sub(s).block(blocknum).b(lead).su(event_plot(j),i).fr_se{types(k)})];
                            patch(patchx(~isnan(patchy)),patchy(~isnan(patchy)),mean([colors(k,:);1,1,1]),'EdgeColor','none','FaceAlpha',0.5)
                            
                            
                            hold on
                            plot(sub(s).block(blocknum).b(lead).su(event_plot(j),i).fr_times{types(k)},sub(s).block(blocknum).b(lead).su(event_plot(j),i).fr{types(k)},'Color',colors(k,:))
                        end
                    end
                    title(event_names(j))
                    xlim([-2,2])
                    ylims=get(gca,'YLim');
                    ylims(1)=0;
                    set(gca,'YLim',ylims);
                end
            end
            suptitle(['Sub',num2str(s),'Block ',num2str(blocknum),'Lead',num2str(lead)])
        end
    end
end
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
%%

types=[1,4,5];
events=[1,2,3,8];
indstruct={};
basestruct={};
for s=1%[1,2,4,7] (output is quite large if you do this for all patients at once, I plan to automate it and do it on cluster eventually)
    for blocknum=1:length(sub(s).block)
        for lead=1:length(sub(s).block(blocknum).b)
            for j=4%length(events)
                for k=1:length(types)
                    if j==4 && k==1
                        continue
                    else
                        trial_type=sub(s).block(blocknum).b(lead).type;
                        trial_type=trial_type==types(k) & cellfun(@(x) ~isempty(x),sub(s).block(blocknum).ao(lead).trials(events(j),:));
                        [wt,f]=cellfun(@(x) cwt(x(7000:13000),1000,'amor','NumOctaves',9),sub(s).block(blocknum).ao(lead).trials(events(j),trial_type),'UniformOutput',0);
                        basemat=cat(3,wt{:});
                        basestruct{s,blocknum,lead,j,k}=basemat;
                    end
                end
            end
        end
    end
end
