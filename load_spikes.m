function spks=load_spikes(file,alg)
spks=struct('t',[],'clust_id',[]);
if strcmp(alg,'wav_clus')
    load(file);
elseif strcmp(alg,'POS')
    num=importdata(file);
    spks.t=num.data(:,3);
    spks.clust_id=num.data(:,2);
    channels=num.data(:,1);
    nC=max(spks.clust_id);
    spks.ch=nan(1,nC);
    for n=1:nC
        spks.ch(n)=mode(channels(spks.clust_id==n));
    end
end
