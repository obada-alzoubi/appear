% automatic cardiac cycle determination directly from EEG-fMRI data
function pkloc = mbp(Fn,srate,acqfrq)
tmbp = tic;

% estimation of cardiac cycle
spt = abs(fft(Fn.*(Fn>=0))).^2;
sgnpt = fft(spt);
tmwin = zeros(size(sgnpt));
hannwin = hann(2*round(8*srate),'periodic');
tmwin(1:round(8*srate)) = hannwin(round(8*srate)+1:2*round(8*srate));
tmwin(length(tmwin)-round(8*srate)+1:length(tmwin)) = hannwin(1:round(8*srate));
sgnpt = sgnpt.*tmwin;
smspt = real(ifft(sgnpt));
smspt = smspt/max(smspt(round(0.5*size(Fn,2)/srate):round(8*size(Fn,2)/srate)));
spt = spt/max(spt(round(0.5*size(Fn,2)/srate):round(8*size(Fn,2)/srate)));
[~, allpklocs] = findpeaks(smspt);
[~, allvalocs] = findpeaks(-smspt);
allvalocs(allvalocs*srate/size(Fn,2)>8) = [];
allpklocs(allpklocs*srate/size(Fn,2)>8) = [];
if(allpklocs(1)<allvalocs(1))
    allpklocs(1) = [];
end
if(allpklocs(end)>allvalocs(end))
    allpklocs(end) = [];
end
allpklocs = zeros(1,numel(allvalocs)-1);
for ii=1:numel(allvalocs)-1
    [~,ind] = max(smspt(allvalocs(ii):allvalocs(ii+1)));
    allpklocs(ii) = ind + allvalocs(ii) - 1;
end

%  find fundamental, bp1, first harmonic and reference frequencies
pkloc = [];
hbr = 120;
while(numel(pkloc)==0 && hbr<=150)
    pkloc = allpklocs(allpklocs>=floor(0.5/srate*length(Fn))&allpklocs<=floor((hbr/60)/srate*length(Fn)));
    hbr = hbr + 7.5;
end
selval = 0;
for jj=1:numel(pkloc)
    valoc = allvalocs(allvalocs<pkloc(jj));
    [~, index] = min(abs(valoc - pkloc(jj)));
    lfvaloc = valoc(index);
    valoc = allvalocs(allvalocs>pkloc(jj));
    [~, index] = min(abs(valoc - pkloc(jj)));
    rgvaloc = valoc(index);
    data = sort(spt(lfvaloc:rgvaloc),'descend');
    tmpv = smspt(pkloc(jj));
    if(tmpv>selval)
        selval = tmpv;
        sptpkloc = pkloc(jj);
        smsplfvaloc = lfvaloc;
        smsprgvaloc = rgvaloc;
    end
end

valoc = allvalocs(allvalocs<=smsplfvaloc-1);
pkloc = allpklocs(allpklocs<=sptpkloc-1);
if(numel(valoc)>0 && numel(pkloc)>0)
    valoc = sort(valoc,'descend');
    pkloc = sort(pkloc,'descend');
    val1 = smspt(sptpkloc);
    val2 = smspt(smsplfvaloc);
    jj = 1;
    while(and(jj<=numel(valoc),jj<=numel(pkloc)) && smspt(pkloc(jj))>=0.15*smspt(sptpkloc) && smspt(pkloc(jj))<val1 && smspt(valoc(jj))<val2 && pkloc(jj)>round(0.5*size(Fn,2)/srate))
        val1 = smspt(pkloc(jj));
        val2 = smspt(valoc(jj));
        smsplfvaloc = valoc(jj);
        jj = jj + 1;
    end
end
rg1 = max(sptpkloc*2-0.5*(smsprgvaloc-smsplfvaloc),smsprgvaloc);
rg2 = min(sptpkloc*2+0.5*(smsprgvaloc-smsplfvaloc),2.5*sptpkloc);
pkloc = allpklocs( allpklocs >= rg1 & allpklocs <= rg2 );
while(numel(pkloc)==0)
    rg1 = max(rg1 - 0.05*(smsprgvaloc-smsplfvaloc),smsprgvaloc);
    rg2 = min(rg2 + 0.05*(smsprgvaloc-smsplfvaloc),2.5*sptpkloc);
    pkloc = allpklocs( allpklocs >= rg1 & allpklocs <= rg2 );
end
[~,ind] = max(smspt(pkloc));
hmpkloc = pkloc(ind);
hmlfvaloc = max(allvalocs(allvalocs<hmpkloc));
valoc = allvalocs(allvalocs>=smsprgvaloc+1 & allvalocs<hmpkloc);
pkloc = allpklocs(allpklocs>=sptpkloc+1 & allpklocs<hmpkloc);
val1 = smspt(sptpkloc);
val2 = smspt(smsprgvaloc);
jj = 1;
while( and(jj<=numel(valoc),jj<=numel(pkloc)) && smspt(pkloc(jj))>=0.15*smspt(sptpkloc) && smspt(pkloc(jj))<val1 && smspt(valoc(jj))<val2 && valoc(jj)<=hmlfvaloc )
    val1 = smspt(pkloc(jj));
    val2 = smspt(valoc(jj));
    smsprgvaloc = valoc(jj);
    rg1 = max(sptpkloc*2-0.5*(smsprgvaloc-smsplfvaloc),smsprgvaloc);
    rg2 = min(sptpkloc*2+0.5*(smsprgvaloc-smsplfvaloc),2.5*sptpkloc);
    pkloc = allpklocs( allpklocs >= rg1 & allpklocs <= rg2 );
    while(numel(pkloc)==0)
        rg1 = max(rg1 - 0.05*(smsprgvaloc-smsplfvaloc),smsprgvaloc);
        rg2 = min(rg2 + 0.05*(smsprgvaloc-smsplfvaloc),2.5*sptpkloc);
        pkloc = allpklocs( allpklocs >= rg1 & allpklocs <= rg2 );
    end
    [~,ind] = max(smspt(pkloc));
    hmpkloc = pkloc(ind);
    hmlfvaloc = max(allvalocs(allvalocs<hmpkloc));
    jj = jj + 1;
end
bp1bpfrq = [smsplfvaloc smsprgvaloc] * srate/size(Fn,2);

lmt = 0.5*size(Fn,2)/srate;
loc1 = sptpkloc;
val = mean(spt(sptpkloc-round(0.05*size(Fn,2)/srate):sptpkloc+round(0.05*size(Fn,2)/srate)));
jj = -1;
rge = [sptpkloc+round((jj*0.05-0.05)*size(Fn,2)/srate),sptpkloc+round((jj*0.05+0.05)*size(Fn,2)/srate)];
while(rge(1)>lmt)
    if(val>mean(spt(rge(1):rge(2)))&&mean(rge)<=smsplfvaloc)
        loc1 = mean(rge);
        val = mean(spt(rge(1):rge(2)));
    end
    jj = jj - 1;
    rge = [sptpkloc+round((jj*0.05-0.05)*size(Fn,2)/srate),sptpkloc+round((jj*0.05+0.05)*size(Fn,2)/srate)];
end
loc1 = min(loc1,smsplfvaloc);
lmt = hmpkloc-round((smsprgvaloc-smsplfvaloc)/2);
loc2 = sptpkloc;
val = mean(spt(sptpkloc-round(0.05*size(Fn,2)/srate):sptpkloc+round(0.05*size(Fn,2)/srate)));
jj = 1;
rge = [sptpkloc+round((jj*0.05-0.05)*size(Fn,2)/srate),sptpkloc+round((jj*0.05+0.05)*size(Fn,2)/srate)];
while(rge(2)<lmt)
    if(val>mean(spt(rge(1):rge(2)))&&mean(rge)>=smsprgvaloc)
        loc2 = mean(rge);
        val = mean(spt(rge(1):rge(2)));
    end
    jj = jj + 1;
    rge = [smsprgvaloc+round((jj*0.05-0.05)*size(Fn,2)/srate),smsprgvaloc+round((jj*0.05+0.05)*size(Fn,2)/srate)];
end
loc2 = max(loc2,smsprgvaloc);
refloc3 = loc2;
refloc1 = loc1;
refloc2 = loc1+0.5*(sptpkloc-loc1);

% calculate bp1
Fnp = Fn;  Fnp(Fnp<0) = 0;  Fnp = Fnp.*abs(Fnp);
ftbp1 = fft(Fnp);
ftbp1(1:floor(bp1bpfrq(1)*size(Fn,2)/srate)) = 0;
ftbp1(length(Fn)-floor(bp1bpfrq(1)*size(Fn,2)/srate)+1:length(Fn)) = 0;
ftbp1(floor(bp1bpfrq(2)*size(Fn,2)/srate)+1:length(Fn)-floor(bp1bpfrq(2)*size(Fn,2)/srate)) = 0;
bp1 = real(ifft(ftbp1));
bp1 = bp1/std(bp1);
[bp1pkloc,bp1lfloc,bp1rgloc,~,~,~,~] = fninfo(bp1);
bp1lfbd = bp1lfloc;
bp1rgbd = bp1rgloc;
nbpk = size(bp1pkloc,2);
ii = 2;
while(1)
    if(ii>=nbpk)
        break;
    end
    prvpk = bp1pkloc(max(ii-10,1):ii-1);
    prvlfloc = bp1lfloc(max(ii-10,1):ii-1);
    prvrgloc = bp1rgloc(max(ii-10,1):ii-1);
    dptmn = mean(bp1(prvpk)-0.5*(bp1(prvlfloc)+bp1(prvrgloc)));
    cond11a = ( bp1(bp1pkloc(ii))-bp1(bp1lfloc(ii))<0.2*dptmn && bp1pkloc(ii)-bp1pkloc(ii-1)<1/(refloc3/size(Fn,2)) );
    cond11b = ( bp1(bp1pkloc(ii))-bp1(bp1rgloc(ii))<0.2*dptmn && bp1pkloc(ii+1)-bp1pkloc(ii)<1/(refloc3/size(Fn,2)) );
    cond11c = ( or(bp1(bp1pkloc(ii))-bp1(bp1lfloc(ii))<0.2*dptmn,bp1(bp1pkloc(ii))-bp1(bp1rgloc(ii))<0.2*dptmn) && bp1rgloc(ii)-bp1lfloc(ii)<1/(refloc3/size(Fn,2)) );
    cond12 = ( bp1pkloc(ii+1)-bp1pkloc(ii-1)<srate/bp1bpfrq(1) );
    if( (cond11a && cond12) || (cond11b && cond12) || (cond11c && cond12) )
        if( cond11a )
            bp1rgbd(ii-1) = bp1pkloc(ii);
            bp1lfbd(ii+1) = bp1pkloc(ii);
            bp1rgloc(ii-1) = bp1rgloc(ii);
        elseif( cond11b )
            bp1lfbd(ii+1) = bp1pkloc(ii);
            bp1rgbd(ii-1) = bp1pkloc(ii);
            bp1lfloc(ii+1) = bp1lfloc(ii);
        elseif( cond11c)
            bp1rgbd(ii-1) = bp1pkloc(ii);
            bp1lfbd(ii+1) = bp1pkloc(ii);
        end
        bp1pkloc(ii) = [];
        bp1lfloc(ii) = [];
        bp1rgloc(ii) = [];
        bp1lfbd(ii) = [];
        bp1rgbd(ii) = [];
        nbpk = nbpk - 1;
        ii = ii - 1;
    end
    ii = ii + 1;
end

%  calculate bp2
[nu,de] = butter(3,[1,acqfrq]/(0.5*srate));
bp2 = filtfilt(nu,de,Fn);
for ii = 2:numel(bp1pkloc)-1
    rge = [bp1lfbd(ii),bp1rgbd(ii)-1];
    bp2(rge(1):rge(2)) = bp2(rge(1):rge(2))/std(bp2(rge(1):rge(2)));
end
rge = [1,bp1lfbd(2)-1];
bp2(rge(1):rge(2)) = bp2(rge(1):rge(2))/std(bp2(rge(1):rge(2)));
rge = [bp1rgbd(numel(bp1pkloc)-1),size(bp2,2)];
bp2(rge(1):rge(2)) = bp2(rge(1):rge(2))/std(bp2(rge(1):rge(2)));
[bp2pkloc,~,~,~,~,~,~] = fninfo(bp2);
bp2pkamp = bp2(bp2pkloc);
bp2lfslp = zeros(size(bp2pkloc));
bp2rgslp = zeros(size(bp2pkloc));
for ii=1:numel(bp2pkloc)
    bp2lfslp(ii) = bp2pkamp(ii)-min(bp2(max(1,bp2pkloc(ii)-round(0.05*srate)):bp2pkloc(ii)));
    bp2rgslp(ii) = bp2pkamp(ii)-min(bp2(bp2pkloc(ii):min(bp2pkloc(ii)+round(0.05*srate),size(bp2,2))));
end
Bp2 = struct('data',bp2,'pkloc',bp2pkloc,'pkamp',bp2pkamp,'lfslp',bp2lfslp,'rgslp',bp2rgslp,'selpkloc',zeros(1,numel(bp1pkloc)));

% select one peak in each estimated cycle
for ii=1:numel(bp1pkloc)
    itvl = [bp1lfbd(ii),bp1rgbd(ii)];
    rge = (max(ii-10,1):ii-1);
    coritvl = zeros(numel(rge),2);
    coritvlpkloc = zeros(1,numel(rge));
    for jj=1:numel(rge)
        coritvl(jj,:) = [bp1lfbd(rge(jj)),bp1rgbd(rge(jj))];
        coritvlpkloc(jj) = Bp2.selpkloc(Bp2.selpkloc>=coritvl(jj,1) & Bp2.selpkloc<coritvl(jj,2));
    end
    [pkind,sel,selflg,Bp2] = selpk2(Bp2,itvl,coritvl,coritvlpkloc,[mean(Bp2.data)+2*std(Bp2.data),mean(max([Bp2.lfslp;Bp2.rgslp],[],1))+2*std(max([Bp2.lfslp;Bp2.rgslp],[],1))],[1 1 1 1]);
    Bp2.selpkloc(ii) = Bp2.pkloc(pkind(find(sel==1)));
end

% add 1st and last peak if needed
nghb = [numel(bp1pkloc),numel(bp1pkloc)-10; 1,10];
bd = [bp1rgbd(end),size(Fn,2);bp1lfbd(1),1];
for ii = 1:size(nghb,1)
    if(max(bd(ii,:))-min(bd(ii,:))>3)
        [~, bp1pkloc1] = findpeaks(bp1(min(bd(ii,:)):max(bd(ii,:))));
        if(numel(bp1pkloc1)>0)
            bp1pkloc1 = bp1pkloc1 + min(bd(ii,:)) - 1;
            [~,ind] = max(bp1(bp1pkloc1));
            bp1pkloc1 = bp1pkloc1(ind);
            Bp1lfrs = bp1(bp1pkloc)-bp1(bp1lfloc);
            Bp1rgrs = bp1(bp1pkloc)-bp1(bp1rgloc);
            dpt = 0.5*(Bp1lfrs(min(nghb(ii,:)):max(nghb(ii,:)))+Bp1rgrs(min(nghb(ii,:)):max(nghb(ii,:))));
            if(1)
                itvl = [min(bd(ii,:)),max(bd(ii,:))];
                rge = nghb(ii,1:2);
                coritvl = zeros(numel(rge),2);
                coritvlpkloc = zeros(1,numel(rge));
                for jj=1:numel(rge)
                    coritvl(jj,:) = [bp1lfbd(rge(jj)),bp1rgbd(rge(jj))];
                    coritvlpkloc(jj) = Bp2.selpkloc(Bp2.selpkloc>=coritvl(jj,1) & Bp2.selpkloc<coritvl(jj,2));
                end
                [pkind,sel,selflg,Bp2] = selpk2(Bp2,itvl,coritvl,coritvlpkloc,[mean(Bp2.data)+2*std(Bp2.data),mean(max([Bp2.lfslp;Bp2.rgslp],[],1))+2*std(max([Bp2.lfslp;Bp2.rgslp],[],1))],[1 1 1 1]);
                if(bd(ii,2)==size(Fn,2) && all(selflg(find(sel==1),:)==1))
                    Bp2.selpkloc = [Bp2.selpkloc,Bp2.pkloc(pkind(find(sel==1)))];
                elseif(bd(ii,2)==1 && all(selflg(find(sel==1),:)==1))
                    Bp2.selpkloc = [Bp2.pkloc(pkind(find(sel==1))),Bp2.selpkloc];
                end
            end
        end
    end
end

% cycle adjustment
sep0 = 0.5*size(Fn,2)/sptpkloc;
sep1 = size(Fn,2)/refloc3;
sep2 = size(Fn,2)/refloc2;
sep3 = size(Fn,2)/refloc1;
ii=2;
flag = 1;
while(flag)
    ii = ii + 1;
    if(ii>=3 && ii<=numel(Bp2.selpkloc)-1)
        dis1 = Bp2.selpkloc(ii)-Bp2.selpkloc(ii-1);
        dis0 = Bp2.selpkloc(ii-1)-Bp2.selpkloc(ii-2);
        dis2 = Bp2.selpkloc(ii+1)-Bp2.selpkloc(ii);
        cond21a = ( dis1<sep1 && or(dis0>sep2,dis2>sep2) && (dis0+dis1+dis2)>=2.5*size(Fn,2)/sptpkloc );
        cond21b = ( dis1<sep1 && ~cond21a );
        cond21c = ( dis1>sep3 );
        if(cond21a || cond21b || cond21c)
            if(cond21c)
                rmind = [];
                isind = ii-1;
                itvl = [Bp2.selpkloc(ii-1)+ceil(sep1),Bp2.selpkloc(ii)];
                if(ii-2>10)
                    rge = (max(ii-1-10,1):ii-2);
                else
                    rge = (1:10);
                end
            elseif(cond21b)
                adind = [];
                itvl = [Bp2.selpkloc(ii-2)+ceil(sep1),Bp2.selpkloc(ii+1)];
                if(ii-2>10)
                    rge = (max(ii-1-10,1):ii-2);
                else
                    rge = (1:10);
                end
            elseif(cond21a)
                [~,ind] = min([sum(abs(diff(Bp2.selpkloc(ii-2:ii)))),sum(abs(diff(Bp2.selpkloc(ii-1:ii+1))))]);
                if(ind==1)
                    rmind = ii;
                    isind = ii;
                    itvl = [Bp2.selpkloc(ii-1)+ceil(sep1),Bp2.selpkloc(ii+1)];
                    if(ii-1>10)
                        rge = (max(ii-10,1):ii-1);
                    else
                        rge = (1:10);
                    end
                elseif(ind==2)
                    rmind = ii-1;
                    isind = ii-1;
                    itvl = [Bp2.selpkloc(ii-2)+ceil(sep1),Bp2.selpkloc(ii-1)];
                    if(ii-2>10)
                        rge = (max(ii-1-10,1):ii-2);
                    else
                        rge = (1:10);
                    end
                end
            end
            coritvl = zeros(numel(rge),2);
            coritvlpkloc = zeros(1,numel(rge));
            for jj=1:numel(rge)
                coritvl(jj,:) = [max(Bp2.selpkloc(rge(jj))-ceil(sep0),1),min(Bp2.selpkloc(rge(jj))+floor(sep0),numel(Bp2.data))];
                coritvlpkloc(jj) = Bp2.selpkloc(rge(jj));
            end
            [pkind,sel,selflg,Bp2] = selpk2(Bp2,itvl,coritvl,coritvlpkloc,[mean(Bp2.data)+2*std(Bp2.data),mean(max([Bp2.lfslp;Bp2.rgslp],[],1))+2*std(max([Bp2.lfslp;Bp2.rgslp],[],1))],[1 1 1 1]);
            if(numel(pkind)>0)
                if(cond21c)
                    adind = pkind(find(sel==1));
                elseif(cond21b)
                    pkrmind = Bp2.pkloc(pkind(find(sel==0)));
                    for jj=[ii-1,ii]
                        if(numel(find(pkrmind==Bp2.selpkloc(ii-1)))>0)
                            rmind = ii-1;
                            isind = rmind;
                        elseif(numel(find(pkrmind==Bp2.selpkloc(ii)))>0)
                            rmind = ii;
                            isind = rmind;
                        else
                            rmind = [];
                            isind = ii;
                        end
                    end
                elseif(cond21a)
                    adind = pkind(find(sel==1));
                end
                Bp2.selpkloc = [Bp2.selpkloc(1:isind),Bp2.pkloc(adind),Bp2.selpkloc(isind+1:end)];
                Bp2.selpkloc(rmind) = [];
                ii = ii - numel(rmind) + numel(adind);
                if(cond21c==1)
                    ii = ii - 1;
                end
            end
        end
    else
        flag = 0;
    end
end

% find corresponding peak in ic
pkloc = zeros(size(Bp2.selpkloc));
for ii=1:numel(Bp2.selpkloc)
    wdt = ceil(srate/(acqfrq-0.5));
    [~,ind] = findpeaks(Fn(max(Bp2.selpkloc(ii)-wdt,1):min(Bp2.selpkloc(ii)+wdt,size(Fn,2))));
    if(numel(ind)>0)
        ind=ind+max(Bp2.selpkloc(ii)-wdt,1)-1;
        [~,ind1] = max(Fn(ind));
        pkloc(ii) = ind(ind1);
    else
        pkloc(ii) = Bp2.selpkloc(ii);
    end
end

fprintf('Multiple scale cardiac peak detection done -- %f second\n',toc(tmbp));
end

