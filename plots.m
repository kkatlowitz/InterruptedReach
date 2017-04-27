colors=get(gca,'colororder');
close gcf

j=3;
for b=1:2
    figure
    
    binl=sub(1).block(b).psth_start_num:0.002:(sub(1).block(b).psth_end_num-sub(1).block(b).binwidth); %sliding window
    binu=(sub(1).block(b).psth_start_num+sub(1).block(b).binwidth):0.002:sub(1).block(b).psth_end_num;
    
    for i=1:size(sub(1).block(b).su,2)
        subplot(2,size(sub(1).block(b).su,2),i)
        plot(sub(1).block(b).su(j,i).stimes_align*1000,sub(1).block(b).su(j,i).strials,'.') %plot rasters
        title(['Single Unit ',num2str(i)])
        xlabel('Time (ms)')
        set(gca,'box','off')
        subplot(2,size(sub(1).block(b).su,2),size(sub(1).block(b).su,2)+i)
        patchx=[(binl+binu)/2*1000,fliplr((binl+binu)/2*1000)];
        patchy=[sub(1).block(b).su(j,i).fr-sub(1).block(b).su(j,i).fr_se,fliplr(sub(1).block(b).su(j,i).fr+sub(1).block(b).su(j,i).fr_se)];
        patch(patchx(~isnan(patchy)),patchy(~isnan(patchy)),mean([colors(1,:);1,1,1]),'EdgeColor','none')
        hold on
        plot((binl+binu)/2*1000,sub(1).block(b).su(j,i).fr)
        title(['Single Unit ',num2str(i)])
        xlabel('Time (ms)')
        set(gca,'box','off')
    end
    suptitle(['Block',num2str(b)])
end

%%
event_plot=[1:6,8];
for b=1:2
    figure
    for j=1:7
        for i=1:size(sub(1).block(b).su,2)
            subplot(7,size(sub(1).block(b).su,2),size(sub(1).block(b).su,2)*(j-1)+i)
            patchx=[(binl+binu)/2*1000,fliplr((binl+binu)/2*1000)];
            patchy=[sub(1).block(b).su(j,i).fr-sub(1).block(b).su(j,i).fr_se,fliplr(sub(1).block(b).su(j,i).fr+sub(1).block(b).su(j,i).fr_se)];
            patch(patchx(~isnan(patchy)),patchy(~isnan(patchy)),mean([colors(1,:);1,1,1]),'EdgeColor','none')
            hold on
            plot((binl+binu)/2*1000,sub(1).block(b).su(j,i).fr)
            line([0,0],[0,100],'color','r')
            if b==1
                ylim([0,100])
            else
                if i==1
                    ylim([0,100])
                else
                    ylim([0,50])
                end
            end
        end
    end
    suptitle(['Block',num2str(b)])
end

%%
event_plot=[1:6,8];
for b=1:2
    figure
    for j=1:7
        for i=1:size(sub(1).block(b).su,2)
            subplot(7,size(sub(1).block(b).su,2),size(sub(1).block(b).su,2)*(j-1)+i)
            [wt,f]=cwt(sub(1).block(b).su(j,i).fr(~isnan(sub(1).block(b).su(j,i).fr)),'morse',500,'WaveletParameters',[3,15],'NumOctaves',8,'Voicesperoctave',20);
            %[wt,f]=cwt(sub(1).block(b).su(j,i).fr(~isnan(sub(1).block(b).su(j,i).fr)),'amor','NumOctaves',8,'Voicesperoctave',20);
            bins=(binu+binl)/2*1000;
            imagesc(abs(wt(51:161,:)))
            xlim([0 1976]);
            set(gca,'XTick',88.5:300:1888.5,'XTickLabel',round(bins(88:300:1888)+1),'YTick',20:20:100,'YTickLabel',round(f(71:20:161)))
        end
    end
end


