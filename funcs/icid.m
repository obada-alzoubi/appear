%This Function is for selecting the artifactual ICs
function [cbicind,saccade_IC,blink_IC,topomap,spectra_BCG_ICs,tpblink,tpsac,smolregion,singchan,muscleloc] = icid(ic,A,mixsig,srte,mriperd)
tol = 1.0e-6;
% Number of ICs is equal to number of channels
szic = size(ic);
szmixmat = size(A);
icnum = szic(1);
szsm = szic(2);
nbch = szmixmat(1);
[nu,de] = butter(3,0.0008,'high');
for ii=1:icnum
    ic(ii,1:szsm) = filtfilt(nu,de,ic(ii,1:szsm));
end
for ii=1:nbch
    mixsig(ii,1:szsm) = filtfilt(nu,de,mixsig(ii,1:szsm));
end

% ic statistics
[icave,icstd,~,~,ickurt] = stat(ic);
regthrs = [0.3,0.1,0.5];
grdn = 100;
thrsnum = size(regthrs,2);
diffr(1) = linintep([5,15],[1/3,0.45],mean(abs(ickurt(1:icnum)))); %Equation S10
diffr(2) = 0.25;
ignsignif = zeros(1,nbch);
ignsignif([17:19,21:24]) = 1;

% spectral analysis
mofrqind = zeros(1,icnum); % motion frequency (CB+RM motions) index
sptcond = zeros(1,icnum);
[spt,frq] = pwelchcal(ic,srte,512,256); %spt is the power spectral density curve, frq is the freq. of the PSD curve
sptsz = size(spt,2);
rmfrqrge = floor([0.5,4.5]/(srte/2)*(sptsz-1))+1; % rapid head motion freq. range
cbfrqrge = floor([2,7]/(srte/2)*(sptsz-1))+1; % cardioballistic freq. range
mofrqrge = floor([0.5,7]/(srte/2)*(sptsz-1))+1; % motion freq. range
nrfrqrge = floor([8,12]/(srte/2)*(sptsz-1))+1; % neuronal freq. range
bkfrqrge = floor([0.5,3]/(srte/2)*(sptsz-1))+1; % blink freq. range
sptscfrqrge = [1,floor(4/(srte/2)*(sptsz-1))+1]; % amplitude scale(S0) freq. range

% Defining different frequency bands
deltaband = floor([0.5,4]/(srte/2)*(sptsz-1))+1;
thetaband = floor([4,8]/(srte/2)*(sptsz-1))+1;
alphaband = floor([8,12]/(srte/2)*(sptsz-1))+1;
betaband = floor([12,30]/(srte/2)*(sptsz-1))+1;
gammaband = floor([30,60]/(srte/2)*(sptsz-1))+1;

% Defining flag for channels with high alpha activity
alphaflg = zeros(1,icnum);
% Defining flag for channels with high alpha activity
gammaflg = zeros(1,icnum);

% Calculating the normalized power ratio different frequency bands

for iIC = 1:icnum
    [wholefrqpk,wholefrqpkloc] = findpeaks(spt(iIC,:));
    [maxpeak,maxpeakloc]=max(wholefrqpk);
    MaxPkLoc=wholefrqpkloc(maxpeakloc);
    
    delta_nm = sum(spt(iIC,deltaband(1):deltaband(2)))/(deltaband(2)-deltaband(1));
    theta_nm = sum(spt(iIC,thetaband(1):thetaband(2)))/(thetaband(2)-thetaband(1));
    alpha_nm = sum(spt(iIC,alphaband(1):alphaband(2)))/(alphaband(2)-alphaband(1));
    beta_nm = sum(spt(iIC,betaband(1):betaband(2)))/(betaband(2)-betaband(1));
    gamma_nm = sum(spt(iIC,gammaband(1):gammaband(2)))/(gammaband(2)-gammaband(1));
    % Select the channels with having higher alpha actvity compared to
    % other bands
    if (alpha_nm > max(max(delta_nm,theta_nm),beta_nm)) || ((alphaband(1) < MaxPkLoc)&& (alphaband(2) > MaxPkLoc))
        alphaflg(1,iIC) = 1;
    else
        alphaflg(1,iIC) = 0;
    end
    % Select the channels with having higher gamma actvity compared to
    % other bands
    if (gamma_nm > max(max(max(delta_nm,theta_nm),beta_nm),alpha_nm))
        gammaflg(1,iIC) = 1; %mark this one as muscle artifact
    end
    
end


sptampsc = zeros(1,icnum);
mofrqpknm = zeros(1,icnum);
mofrqpkrs = zeros(1,icnum);
cbfrqpkrs = zeros(1,icnum);
maxcbfrqpk = zeros(1,icnum);
avepowrm = zeros(1,icnum);
avepowcb = zeros(1,icnum);
rmconvwdt = zeros(1,icnum);
rmconvdpt = zeros(1,icnum);
mrfrqpkrs = zeros(1,icnum);
lgmopk = zeros(1,icnum);
condbkspt = zeros(1,icnum);
condscspt = zeros(1,icnum);
condlgmopkspt = zeros(1,icnum);
condrspspt = zeros(1,icnum);

% spectral peaks identification
% S0: sptampsc
% rnr: mrfrqpkrs
% pkrse: rcb the peak rise in the CB frequency range
for ii = 1:icnum
    [mofrqpk,mofrqpkloc] = findpeaks(spt(ii,mofrqrge(1):mofrqrge(2))); % find the spectrum peaks in MO freq. range in the PSD curve
    mofrqpknm(ii) = numel(mofrqpk); % how many MO peaks are there in each IC
    if(mofrqpknm(ii)>0) % If there is a peak...
        mofrqpkloc = mofrqpkloc + mofrqrge(1) - 1; % The MO peak loc + lower bound of MO freq. range - 1 (correct the peak location)
    end
    [nrfrqpk,nrfrqpkloc] = max(spt(ii,nrfrqrge(1):nrfrqrge(2))); % find the spectrum max in neuronal freq. range in the PSD curve
    nrfrqpkloc = nrfrqpkloc + nrfrqrge(1) - 1; % The NR peak loc + lower bound of NR freq. range - 1 (correct the peak location)
    [pkloc,~,~,~,dpt,~,~] = fninfo(spt(ii,1:sptsz)); % 0.05*S0
    sptampsc(ii) = max(spt(ii,sptscfrqrge(1):sptscfrqrge(2)))-min(spt(ii,sptscfrqrge(1):sptscfrqrge(2))); % S0
    minpkamp = 0.05*sptampsc(ii); % 0.05*S0
    dscdpk1 = zeros(1,mofrqpknm(ii)); % discard peak locations; zeros the size of number of peaks found in the MO freq. range
    for jj=1:mofrqpknm(ii) %1 to the number of peaks
        dpt1 = dpt(pkloc==mofrqpkloc(jj)); %where Vr is below 8Hz
        if(dpt1<minpkamp) %if the value is less than 0.05*S0...
            dscdpk1(jj) = 1; %mark the location to be discarded
        end
    end
    mofrqpknm(ii) = mofrqpknm(ii) - sum(dscdpk1); %update the number of MO peaks after some have been discarded
    mofrqpkloc = mofrqpkloc(dscdpk1==0); %update the locations of MO peaks after some have been discarded
    
    [~,minimaloc] = findpeaks(-spt(ii,1:sptsz)); %find location of all local minima of the PSD curve
    [nrpklftmin,~] = min(abs(minimaloc(minimaloc < nrfrqrge(1))-nrfrqrge(1))); % The difference between the start of NR range and the last minimal before starting the NR range
    %Table S1 aEREMCOR Peak rise in NR range
    if(numel(nrpklftmin)>0) %if the NR peak minimum exists then..
        nrpklftmin = nrfrqrge(1) - nrpklftmin; %loc of the left minima from NR
        mrfrqpkrs(ii) = nrfrqpk - min(spt(ii,nrpklftmin:nrfrqpkloc)); %if Vl exists
    elseif(numel(nrpklftmin)==0)%if Vl doesn't exist
        mrfrqpkrs(ii) = nrfrqpk - min(spt(ii,nrfrqrge(1):nrfrqpkloc));
    end
    maxcbfrqpk(ii) = min(spt(ii,cbfrqrge(1):cbfrqrge(2))); %minimum value found in the CB freq. range (in cb range, location needs to be corrected since it is before shifting)
    nbkpkrs = 0;
    nbkpkflg = 0;
    cbpkflg = 0;
    rgtminflg = 0;
    for jj=1:mofrqpknm(ii) %for 1:number of MO peaks
        [rgtmin,~] = min(abs(minimaloc(minimaloc > mofrqpkloc(jj)) - mofrqpkloc(jj))); %right min; smallest distance above the local minima from the MO peak
        [lftmin,~] = min(abs(minimaloc(minimaloc < mofrqpkloc(jj)) - mofrqpkloc(jj))); %left min; smallest distance below the local minima from the MO peak
        if( numel(lftmin)==0 ) %if no left min
            lftmin = 1; %make left min = 1
            lftminflg = 0;
        else
            lftmin = mofrqpkloc(jj) - lftmin; %Adjust the correct loc. left min is the loc, after removing left min from the MO peak
            lftminflg = 1;
        end
        rgtmin = mofrqpkloc(jj) + rgtmin; %Adjust the correct loc.right min is the loc after adding right min to the MO peak
        if( rgtminflg==0 && frq(rgtmin)<=7 ) %if local minimum is below 7 Hz
            rgtminflg = 1;
        else
            rgtminflg = 0;
        end
        
        lftrse = spt(ii,mofrqpkloc(jj))-spt(ii,lftmin); %left rise; PSD value of MO peak - left local minima
        rgtrse = spt(ii,mofrqpkloc(jj))-spt(ii,rgtmin); %right rise; PSD value of MO peak - right local minima
        pkrse = 0.5*(lftrse+rgtrse); %peak rise
        if( lftminflg==1 && frq(mofrqpkloc(jj))>0.5 && frq(mofrqpkloc(jj))<4.5 ) % If the motion peak is in the RM range
            if( min(lftrse,rgtrse)>minpkamp )
                nbkpkrs = max( nbkpkrs, pkrse );
                nbkpkflg = 1;
            end
        end
        if( lftminflg==1 && frq(mofrqpkloc(jj))>2 && frq(mofrqpkloc(jj))<7 ) % If the motion peak is in the CB range
            if( pkrse>0.2*sptampsc(ii) ) % rcb>0.2*S0 condition (3) Large MO peak
                cbfrqpkrs(ii) = max(cbfrqpkrs(ii),pkrse);% max(rcb)>0.2S0
                maxcbfrqpk(ii) = max(maxcbfrqpk(ii),spt(ii,mofrqpkloc(jj))); % maxcbfrqpk=max (peak cb & mo)
                cbpkflg = 1;
            end
        end
        mofrqpkrs(ii) = max(mofrqpkrs(ii),pkrse); %max peak rise in cb range
    end
    avepowrm(ii) = mean(spt(ii,rmfrqrge(1):rmfrqrge(2)))-min(spt(ii,1:nrfrqpkloc)); % mean(Srm)-Smin
    avepowcb(ii) = mean(spt(ii,cbfrqrge(1):cbfrqrge(2)))-min(spt(ii,1:nrfrqpkloc)); % mean(Scb)-Smin
    
    % spectrum slope in motion freq range
    rmsptsc = max(spt(ii,rmfrqrge(1):rmfrqrge(2)))-min(spt(ii,rmfrqrge(1):rmfrqrge(2)));
    rmsptnrm = ( spt(ii,1:sptsz) - min(spt(ii,rmfrqrge(1):rmfrqrge(2))) ) / rmsptsc; % Normalize the spectrum in RM freq. range
    deri1 = diff(rmsptnrm(1:sptsz)); % the first derivative of normalized spectrum
    smslp = deri1;
    % smooth the local fluctuation of each point by taking the average
    % with its neighbors
    for jj=2:sptsz-2
        smslp(jj) = mean(deri1(jj-1:jj+1));
    end
    deri2 = diff(smslp); % the second derivative of normalized spectrum
    trnptflg = zeros(size(deri2));
    for jj=2:sptsz-3
        if(deri2(jj-1)>=0 && deri2(jj+1)<=0 && deri2(jj-1)-deri2(jj+1)>0.02) % S"(v0-)>=0, S"(V0+)<=0 and S"(V0-)-S"(V0+)>0.02
            trnptflg(jj) = 1;
        end
    end
    trnpt = find(trnptflg==1); % Find reflection points
    % Find reflection points in RM range
    trnpt = trnpt(trnpt>rmfrqrge(1));
    trnpt = trnpt(trnpt<rmfrqrge(2));
    for jj=1:numel(trnpt)
        if(jj<numel(trnpt))
            rgtpt = trnpt(jj+1);
        elseif(jj==numel(trnpt)) %last reflection point
            [~,rgtpt] = min(rmsptnrm(trnpt(jj):rmfrqrge(2)));
            rgtpt = rgtpt + trnpt(jj) - 1;
        end
        tmpwdt = (rgtpt-trnpt(jj))*(srte/2)/(sptsz-1);
        tmpdpt = rmsptnrm(trnpt(jj))-rmsptnrm(rgtpt); %Rrm'=S(V0)-S(V'0)
        if( tmpdpt>rmconvdpt(ii) )
            rmconvwdt(ii) = tmpwdt;
            rmconvdpt(ii) = tmpdpt;
        end
    end
    
    % spectrum slope in blink freq range
    bksptsc = max(spt(ii,bkfrqrge(1):bkfrqrge(2))) - min(spt(ii,bkfrqrge(1):bkfrqrge(2)));
    bksptnrm = ( spt(ii,bkfrqrge(1):bkfrqrge(2)) - min(spt(ii,bkfrqrge(1):bkfrqrge(2))) ) / bksptsc;
    bkflt = std(diff(diff(bksptnrm)));
    nbkpkind = 0;
    if(nbkpkflg==1 && nbkpkrs>2/3*sptampsc(ii))
        nbkpkind = 1;
    end
    condbkspt(ii) = or( bkflt<0.2 && not(nbkpkind), mofrqpknm(ii)==0 );
    condscspt(ii) = condbkspt(ii);
    
    % selection
    sussptcond = 0;
    if(mofrqpknm(ii)>=1)
        sussptcond = ( lftminflg==0 && rgtminflg==0 );
    end
    minrmfrqval = min(spt(ii,rmfrqrge(1):rmfrqrge(2)));
    avemopkval = mean(spt(ii,mofrqpkloc));
    lgmopk(ii) = ( mofrqpkrs(ii)>0.2*sptampsc(ii) );
    smmopk = ( not(lgmopk(ii)) && avemopkval>nrfrqpk-2 ); %Condition (4) Small MO peaks
    smnrpk = or( mrfrqpkrs(ii)<0.2*sptampsc(ii), minrmfrqval>nrfrqpk-2 ); %Condition (5) Small NR peak
    domnrpk = ( mrfrqpkrs(ii)>sptampsc(ii)/3 ); % If rise in NR range is larger than 0.33 of the largest rise in all spectrum range
    if(domnrpk==1)
        cond1 = ( cbpkflg==1 && cbfrqpkrs(ii)>mrfrqpkrs(ii)-3 ); %Condition (iii):max(rcb)>rnr-3
        cond2 = ( cbpkflg==1 && maxcbfrqpk(ii)>nrfrqpk-3 && avepowcb(ii)>mrfrqpkrs(ii)/3 ); % Condition (iv)maxcbfrqpk: max(Scb) and nrfrqpk: max(Snr)
        cond3 = ( avepowrm(ii)>1.75*mrfrqpkrs(ii) ); %condition (ii) in aEremcor mean(S(v)-Smin>a2rnr
        accnrpk = ( cond1 || cond2 || cond3 );
    else
        accnrpk = 1;
    end
    if( accnrpk && not(sussptcond) )
        if( mofrqpknm(ii)>=1 && lgmopk(ii) )
            mofrqind(ii) = 1;
            if( cbpkflg==1 )
                condlgmopkspt(ii) = 1;
            end
        elseif( mofrqpknm(ii)>=1 && and(smmopk,smnrpk) )
            mofrqind(ii) = 1;
        elseif( rmconvdpt(ii)>0.15 && smnrpk )  %Condition 6: obvious negative convexity
            mofrqind(ii) = 1;
        end
    end
    condrspspt(ii) = accnrpk;
end

% topo map analysis
motpind = zeros(1,icnum);
polcond = zeros(1,icnum);
cbtpind = zeros(1,icnum);
[xq, yq] = meshgrid(-1:2/grdn:1, -1:2/grdn:1);
regnm = zeros(icnum,thrsnum,2);
regArcnm = zeros(icnum,thrsnum,2);
regA = zeros(icnum,thrsnum,2);
regArcA = zeros(icnum,thrsnum,2);
regcntval = zeros(icnum,thrsnum,2);
regcnt2org = zeros(icnum,thrsnum);
reg12cnt2org = zeros(icnum,1);
regX = zeros(icnum,thrsnum,2);
regY = zeros(icnum,thrsnum,2);
regR = zeros(icnum,thrsnum,2);
regT = zeros(icnum,thrsnum,2);
rregnm = zeros(icnum,thrsnum);
rregArcnm = zeros(icnum,thrsnum);
rregA = zeros(icnum,thrsnum);
rregArcA = zeros(icnum,thrsnum);
condbktp = zeros(1,icnum);
condsctp = zeros(1,icnum);
condrmmot = zeros(1,icnum);
sustpcond = zeros(1,icnum);
x = zeros(1,nbch);
y = zeros(size(x));
% Channel Location for 32-channel BP EEG system
chpos = [[-90,-72];[90,72];[-60,-51];[60,51];[-45,0];[45,0];[-60,51];[60,-51];...
    [-90,72];[90,-72];[-90,-36];[90,36];[-90,0];[90,0];[-90,36];[90,-36];...
    [45,90];[0,0];[45,-90];[90,-90];[-31,-46];[31,46];[-31,46];[31,-46];...
    [-69,-21];[69,21];[-69,21];[69,-21];[-113,18];[113,-18];[67,-90]];

for ii=1:nbch
    x(ii) = chpos(ii,1)*pi/180 * cos(chpos(ii,2)*pi/180);
    y(ii) = chpos(ii,1)*pi/180 * sin(chpos(ii,2)*pi/180);
end
x1 = x./max(sqrt(x.^2+y.^2));
y1 = y./max(sqrt(x.^2+y.^2));
x = x1;
y = y1;
toosmall = zeros(2,size(ic,1)); % For the polarity, define the small regions to be removed later
pertop = zeros(size(ic,1),3);
%   perbot = zeros(size(ic,1),3);
for ii=1:icnum
    v(1:nbch) = A(1:nbch,ii) / max(abs(A(1:nbch,ii)));
    vq = griddata(x,y,v,xq,yq,'v4');
    vq = vq / max(max(abs(vq)));
    %% added
    if ii==2;
        ahmad=1;
    end
    vq(xq.^2 + yq.^2 > 1) = 0; %Topographic map in a circle imagesc(vq)
    vqarc = vq;
    %imagesc(vq)
    vqarc(sqrt(xq.^2 + yq.^2) < 0.8) = 0; %Map boundry of 0.2 (arc ciecle)
    if(ii==1)
        bdmap = zeros(size(vq));
        bdinfo = bwboundaries(vq,8,'noholes');
        for ll=1:length(bdinfo{1})
            bdmap(bdinfo{1}(ll,1),bdinfo{1}(ll,2)) = 1; %bdmap is the outer circle of the map
        end
    end
    totA = numel(vq(xq.^2+yq.^2<=1)); %The number of pixels in outer circle
    totArcA = numel(vqarc(vqarc~=0)); %The number of pixels in the arc (the area between outer circle and inner circle)
    
    % regions definition
    for jj=1:thrsnum % jj = 1,2 are thresholds for the polarity regions (page 3 of supplementary for aE-REMOCOR); jj=1 (K=0.3) are primary polarity regions, jj=2 (K'=0.1) are secondary polarity regions
        map = ones(size(vq));
        
        %added by Kaylee
        squaretopmap = zeros(size(vq));
        squarebotmap = zeros(size(vq));
        emptymap = ones(size(vq));
        emptymap(xq.^2 + yq.^2 > 1) = 0;
        squaretopmap(1:20,20:80)=1;
        squarebotmap(80:100,20:80)=1;
        mapbaby_top = emptymap.*squaretopmap;
        mapbaby_bot = emptymap.*squarebotmap;
        mapbaby = mapbaby_top + mapbaby_bot;
        
        map(xq.^2 + yq.^2 > 1) = 0;
        map(vq>regthrs(jj))=0; %positive polarity
        map(vq<-regthrs(jj))=0; %negative polarity
        maplabl = bwlabel(map,8); %separate founded regions (neutral ones |g(r)|<thrdnum(K)
        selregR = sign(maplabl);
        rregnm(ii,jj) = max(max(maplabl)); % Count number of brain regions for each IC(ii) and each thrsnum(jj); column 1 is # of primary polarity regions, column 2 is # of secondary polarity regions, column 3 is # of other regions?
        npatA = zeros(1,rregnm(ii,jj));
        dscdR = 0;
        ignR = 0;
        for ll=1:rregnm(ii,jj)
            npatA(ll) = sum(sign(maplabl(maplabl==ll)))/totA; %The percentage of neutral area to the area of outer circle
            if(npatA(ll)<=0.04 && npatA(ll)>0.001 && dscdR<2) % If the region is very small, just discard it
                dscdR = dscdR + 1; % Add it to Discarded neutral  brain regions
                selregR(maplabl==ll) = 0;
            end
            if(npatA(ll)<=0.001) % If the region is very small, just discard it
                ignR = ignR + 1;
                selregR(maplabl==ll) = 0;
            end
        end
        rregnm(ii,jj) = rregnm(ii,jj) - dscdR - ignR; %number of neutral regions after discarding/ignoring
        map(selregR==0) = 0; %removes the discard/ignore neutral regions from the map
        arcmap = map;
        arcmap(vqarc==0) = 0; %overlap between neutral region map and the arc
        rregArcnm(ii,jj) = max(max(bwlabel(arcmap,8)));
        rregA(ii,jj) = numel(map(map~=0))/totA;  %The percentage of neutral area to the area of outer circle
        rregArcA(ii,jj) = numel(arcmap(arcmap~=0))/totArcA;  %The percentage of neutral area in the arc to the area of outer circle
        for kk=1:2 %1 is positive and 2 negative regions
            ignmap = zeros(size(vq));
            ignmap((-1)^(kk+1)*vq<=regthrs(jj)) = 1;
            ignmap(xq.^2 + yq.^2 > 1) = 0;
            ignmaplabl = bwlabel(ignmap,8);
            ignpatA = zeros(1,max(max(ignmaplabl)));
            for ll=1:max(max(ignmaplabl))
                ignpatA(ll) = sum(sign(ignmaplabl(ignmaplabl==ll)))/totA; %percentage ratio of anything outside the ll region
                %discarding small regions
                if(ignpatA(ll)<=0.025 && ignpatA(ll)>0.001)
                    ignmaplabl(ignmaplabl==ll) = 0;
                end
                if(ignpatA(ll)<=0.001)
                    ignmaplabl(ignmaplabl==ll) = 0;
                end
            end
            map = vq;
            map((-1)^(kk+1)*vq>regthrs(jj)) = 1;
            map((-1)^(kk+1)*vq<=regthrs(jj)) = 0;
            maplabl = bwlabel(map,8);
            selreg = sign(maplabl);  %selected region (one of the positive or negative region)
            
            
            together=mapbaby.*selreg;
            pertop(ii,jj,kk) = sum(sum(together))/sum(sum(mapbaby));

            regnm(ii,jj,kk) = max(max(maplabl));
            patA = zeros(1,regnm(ii,jj,kk));
            dscd = 0;
            ign = 0;
            arcmap = sign(abs(vqarc)); %part of this region that is inside the arc
            arcmap(ignmaplabl>0) = 0;
            
            for ll=1:regnm(ii,jj,kk)
                patA(ll) = sum(sign(maplabl(maplabl==ll)))/totA; % The percentage ratio of the ll region over the area f outer circle
                %If it's too small, ignore it
                if(patA(ll)<=0.025 && patA(ll)>0.001 && dscd<2)
                    if kk==1
                        toosmall (1,ii) = toosmall(1,ii)+patA(ll);
                    else
                        toosmall (2,ii) = toosmall(2,ii)+patA(ll);
                    end
                    dscd = dscd + 1;
                    selreg(maplabl==ll) = 0;
                    arcmap(maplabl==ll) = 0;
                end
                if(patA(ll)<=0.001)
                    ign = ign + 1;
                    selreg(maplabl==ll) = 0;
                    arcmap(maplabl==ll) = 0;
                end
            end
            regnm(ii,jj,kk) = regnm(ii,jj,kk) - dscd - ign; %number of neutral regions after discarding/ignoring
            arcmap(selreg==0) = 0;
            regArcnm(ii,jj,kk) = max(max(bwlabel(arcmap,8)));
            
            regA(ii,jj,kk) = numel(vq(selreg>0))/totA; %percentage of the positive/negative region area to the entire circle area
            regArcA(ii,jj,kk) = numel(vqarc((-1)^(kk+1)*vqarc>regthrs(jj)&selreg>0))/totArcA; %percentage of positive/negative region arc area to the entire arc area
            regX(ii,jj,kk) = mean(xq(selreg>0));
            regY(ii,jj,kk) = mean(yq(selreg>0));
            regR(ii,jj,kk) = mean(sqrt((xq(selreg>0)).^2+(yq(selreg>0)).^2)); % The radius of the center of positive/negative region
            regT(ii,jj,kk) = atan(abs(regY(ii,jj,kk))/abs(regX(ii,jj,kk)));
            % Fixing the Theta
            if(regX(ii,jj,kk)<0 && regY(ii,jj,kk)>0)
                regT(ii,jj,kk) = pi - regT(ii,jj,kk);
            elseif(regX(ii,jj,kk)<0 && regY(ii,jj,kk)<0)
                regT(ii,jj,kk) = regT(ii,jj,kk) + pi;
            elseif(regX(ii,jj,kk)>0 && regY(ii,jj,kk)<0)
                regT(ii,jj,kk) = 2*pi - regT(ii,jj,kk);
            end
            if(regArcA(ii,jj,kk)>tol)
                regX(ii,jj,kk) = regR(ii,jj,kk)*cos(regT(ii,jj,kk));
                regY(ii,jj,kk) = regR(ii,jj,kk)*sin(regT(ii,jj,kk));
            end
            [~, locx] = min(abs((-1:2/grdn:1)-regX(ii,jj,kk)));
            [~, locy] = min(abs((-1:2/grdn:1)-regY(ii,jj,kk)));
            regcntval(ii,jj,kk) = vq((locx-1)*(grdn+1)+locy); %find the gradient value at the x,y location
        end
        regcnt2org(ii,jj) = lne2pt([regX(ii,jj,1) regY(ii,jj,1) 0],[regX(ii,jj,2) regY(ii,jj,2) 0],[0 0 0]); %finds distance d
    end
    
    % blink identification
    bkthrsapp = 1;
    condbktp1=0;
    condbktp2=0;
    condbktp3=0;
    [~, nseposx] = min(abs((-1:2/grdn:1)-0));
    [~, nseposy] = min(abs((-1:2/grdn:1)-1));
    nseposval = vq((nseposx-1)*(grdn+1)+nseposy);
    for kk=1:2
        if( (-1)^(kk+1)*nseposval>2/3 )
            bkregsgn = kk;
            map = vq;
            map((-1)^(bkregsgn+1)*vq>regthrs(bkthrsapp)) = 1;
            map((-1)^(bkregsgn+1)*vq<=regthrs(bkthrsapp)) = 0;
            mapnse = bwselect(map,nseposx,nseposy,8);
            arcmapnse = mapnse;
            arcmapnse(vqarc==0) = 0;
            arcmapnsenm = max(max(bwlabel(arcmapnse,8)));
            smapnse = mapnse;
            smapnse(abs(vq)<2/3) = 0;
            smapnsenm = max(max(bwlabel(smapnse,8)));
            sArcmapnse = smapnse;
            sArcmapnse(vqarc==0) = 0;
            bkArcA = sum(sum(arcmapnse==1))/totArcA;
            bkA = sum(sum(mapnse==1))/totA;
            bksArcA = sum(sum(sArcmapnse==1))/totArcA;
            bksA = sum(sum(smapnse==1))/totA;
            bkbdinfo = bwboundaries(mapnse,8,'noholes');
            bkbd = max(length(bkbdinfo{1})/floor(grdn/2)-2*pi*bkArcA,0);
            rge = [min(((-1)^bkregsgn)*[-0.3,0.45]),max(((-1)^bkregsgn)*[-0.3,0.45])];
            nbkA = numel(vq(xq.^2+yq.^2<=1 & vq>=rge(1) & vq<=rge(2)))/totA;
            bkregX = mean(xq(mapnse==1));
            bkregY = mean(yq(mapnse==1));
            bkregR = mean(sqrt((xq(mapnse==1)).^2+(yq(mapnse==1)).^2));
            bkregT = atan(abs(bkregY)/abs(bkregX));
            if(bkregX<0 && bkregY>0)
                bkregT = pi - bkregT;
            elseif(bkregX<0 && bkregY<0)
                bkregT = bkregT + pi;
            elseif(bkregX>0 && bkregY<0)
                bkregT = 2*pi - bkregT;
            end
            if(bkArcA>tol)
                bkregX = bkregR*cos(bkregT);
                bkregY = bkregR*sin(bkregT);
            end
            
            Anglim = pi/2 - 4/9*pi;
            [~, loc1x] = min(abs((-1:2/grdn:1)-cos(Anglim)));
            [~, loc1y] = min(abs((-1:2/grdn:1)-sin(Anglim)));
            [~, loc2x] = min(abs((-1:2/grdn:1)+cos(Anglim)));
            [~, loc2y] = min(abs((-1:2/grdn:1)-sin(Anglim)));
            condbktp1 = ( arcmapnsenm==1 && smapnsenm==1 && bkregY>0 && abs(bkregX)<bkregY*tan(35*pi/180) && smapnse(loc1y,loc1x)==0 && smapnse(loc2y,loc2x)==0 );
            if(condbktp1)
                secAlim = zeros(1,2);
                condbktp2a = ( bkArcA>13/72 && bkArcA<4/9 && bksArcA>0.45*bkArcA );
                SectorArea = ( 2*pi*bkArcA-sin(2*pi*bkArcA) )/2 /pi;
                secAlim(1:2) = [2,0.9].*( 2*pi*[13/72,4/9]-sin(2*pi*[13/72,4/9]) )/2 /pi;
                Alim1 = max(0.8*SectorArea, secAlim(1));
                Alim2 = min(SectorArea+0.8*sin(bkArcA*2*pi)/2/pi,secAlim(2));
                condbktp2b = ( bkA>Alim1 && bkA<Alim2 && bksA>0.25*bkA );
                sectbd = 1.75 * 2*sin(pi*bkArcA);
                condbdtp2c = ( bkbd<sectbd );
                contbdtp2d = ( nbkA>0.6 || and(nbkA<0.6,not(condlgmopkspt(ii))) );
                condbktp2 = ( condbktp2a && condbktp2b && condbdtp2c && contbdtp2d );
            end
            [~, locx] = min(abs((-1:2/grdn:1)-regX(ii,bkthrsapp,kk)));
            [~, locy] = min(abs((-1:2/grdn:1)-regY(ii,bkthrsapp,kk)));
            condbktp3 = ( mapnse(locy,locx)==1 );
        end
    end
    if( condbktp1 && condbktp2 && condbktp3 )
        condbktp(ii) = 1;
    end
    if( rregnm(ii,3)<=1 && all(eq(regArcnm(ii,3,1:2),1)) && all(eq(regnm(ii,3,1:2),1)) )
        [~,minind] = min(regA(ii,3,1:2));
        maxind = find([1,2]~=minind);
        sustpcond(ii) = ( regY(ii,3,minind)>0 && regY(ii,3,maxind)<0 && ...
            max(regA(ii,3,1:2))>0.45 && min(regA(ii,3,1:2))<0.25 && ...
            abs(regX(ii,3,minind))<abs(regY(ii,3,minind))*tan(35*pi/180) );
    end
    
    % saccade identification
    scthrsapp = 1;
    condsctp1=0;
    condsctp2=0;
    condsctp3=0;
    [~, scpos1x] = min(abs((-1:2/grdn:1)-(-0.9*cos(35/180*pi))));
    [~, scpos1y] = min(abs((-1:2/grdn:1)-0.9*sin(35/180*pi)));
    scpos1val = vq((scpos1x-1)*(grdn+1)+scpos1y);
    [~, scpos2x] = min(abs((-1:2/grdn:1)-0.9*cos(35/180*pi)));
    [~, scpos2y] = min(abs((-1:2/grdn:1)-0.9*sin(35/180*pi)));
    scpos2val = vq((scpos2x-1)*(grdn+1)+scpos2y);
    scArcA = zeros(1,2);
    scA = zeros(1,2);
    scsArcA = zeros(1,2);
    scsA = zeros(1,2);
    scbd = zeros(1,2);
    scregsgn = zeros(1,2);
    scregX = zeros(1,2);
    scregY = zeros(1,2);
    scregR = zeros(1,2);
    scregT = zeros(1,2);
    for kk=1:2
        if( (-1)^(kk+1)*scpos1val>regthrs(scthrsapp) && (-1)^(kk+1)*scpos2val<-regthrs(scthrsapp) )
            scregsgn(1) = kk;
            scregsgn(2) = find([1,2]~=scregsgn(1));
            map = vq;
            map((-1)^(scregsgn(1)+1)*vq>regthrs(scthrsapp)) = 1;
            map((-1)^(scregsgn(1)+1)*vq<=regthrs(scthrsapp)) = 0;
            mapsc1 = bwselect(map,scpos1x,scpos1y,8);
            map((-1)^(scregsgn(2)+1)*vq>regthrs(scthrsapp)) = 1;
            map((-1)^(scregsgn(2)+1)*vq<=regthrs(scthrsapp)) = 0;
            mapsc2 = bwselect(map,scpos2x,scpos2y,8);
            nscA = numel(vq(xq.^2+yq.^2<=1 & vq>=-0.3 & vq<=0.3))/totA;
            for kkp=1:2
                if(kkp==1)
                    mapsc = mapsc1;
                elseif(kkp==2)
                    mapsc = mapsc2;
                end
                scArcA(kkp) = numel(vqarc(vqarc~=0 & mapsc==1))/totArcA;
                scA(kkp) = sum(sign(abs(vq(mapsc==1))))/totA;
                scsArcA(kkp) = sum(sign(abs(vqarc(vqarc~=0 & mapsc==1 & abs(vq)>0.5))))/totArcA;
                scsA(kkp) = sum(sign(abs(vq(mapsc==1 & abs(vq)>0.5))))/totA;
                scbdinfo = bwboundaries(mapsc,8,'noholes');
                scbd(kkp) = max(length(scbdinfo{1})/floor(grdn/2)-2*pi*scArcA(kkp),0);
                if(scsA(kkp)>0)
                    scregX(kkp) = mean(xq(mapsc==1 & abs(vq)>0.5));
                    scregY(kkp) = mean(yq(mapsc==1 & abs(vq)>0.5));
                else
                    scregX(kkp) = mean(xq(mapsc==1));
                    scregY(kkp) = mean(yq(mapsc==1));
                end
                scregR(kkp) = mean(sqrt((xq(mapsc==1)).^2+(yq(mapsc==1)).^2));
                scregT(kkp) = atan(abs(scregY(kkp))/abs(scregX(kkp)));
                if(scregX(kkp)<0 && scregY(kkp)>0)
                    scregT(kkp) = pi - scregT(kkp);
                elseif(scregX(kkp)<0 && scregY(kkp)<0)
                    scregT(kkp) = scregT(kkp) + pi;
                elseif(scregX(kkp)>0 && scregY(kkp)<0)
                    scregT(kkp) = 2*pi - scregT(kkp);
                end
                if(scArcA(kkp)>tol)
                    scregX(kkp) = scregR(kkp)*cos(scregT(kkp));
                    scregY(kkp) = scregR(kkp)*sin(scregT(kkp));
                end
            end
            
            dis = lne2pt([scregX(1),scregY(1),0],[scregX(2),scregY(2),0],[0,0,0]);
            incnt = atan(abs(scregY(2)-scregY(1))/abs(scregX(2)-scregX(1)));
            condsctp1 = ( dis>0.3 && incnt<35/180*pi );
            if(condsctp1)
                condsctp2a = ( min(scArcA)<0.3 && max(scArcA)<0.4 && sum(scArcA)<2/3 );
                condsctp2b = ( max(scA)<0.4 && max(scsA./scA)>0.3 && all(gt(scsArcA./scArcA,0.9*scsA./scA)+eq(scsArcA./scArcA,0.9*scsA./scA)) );
                sectbd = 4*sin(pi*scArcA);
                condsctp2c = ( scbd(1)<sectbd(1) && scbd(2)<sectbd(2) );
                condsctp2d = ( nscA>0.45 );
                condsctp2e = ( abs(scA(1)-scA(2))/max(scA)<2/3 );
                condsctp2 = ( condsctp2a && condsctp2b && condsctp2c && condsctp2d && condsctp2e);
            end
            [~, locx] = min(abs((-1:2/grdn:1)-regX(ii,scthrsapp,scregsgn(1))));
            [~, locy] = min(abs((-1:2/grdn:1)-regY(ii,scthrsapp,scregsgn(1))));
            condsctp3a = ( mapsc1(locy,locx)==1 );
            [~, locx] = min(abs((-1:2/grdn:1)-regX(ii,scthrsapp,scregsgn(2))));
            [~, locy] = min(abs((-1:2/grdn:1)-regY(ii,scthrsapp,scregsgn(2))));
            condsctp3b = ( mapsc2(locy,locx)==1 );
            condsctp3 = ( condsctp3a && condsctp3b );
        end
    end
    if( condsctp1 && condsctp2 && condsctp3 )
        condsctp(ii) = 1;
    end
    
    % selection
    % kaylee added condition 5:
    cond5 = 0;
    for kk=1:2
        kkp = find([1,2]~=kk);
        if (pertop(ii,1,kk)>=0.41 && pertop(ii,1,kkp)>=0.41)
            cond5 = 1;
            break;
        end
    end
    
    cond1 = ( rregnm(ii,1)<=1 ); %Condition 1: no more than one neutral region in the topographic map
    % cond2a: only one positive/negative arc and and region
    cond2a = ( all(eq(regArcnm(ii,1,1:2),1)) && all(eq(regnm(ii,1,1:2),1)) ); %Condition 2: Only one positive/negative polarity region and polarity arc region is allowed
    cond2b = 0;
    for kk=1:2 %for two polarity regions
        kkp = find([1,2]~=kk);
        if( regArcnm(ii,1,kk)==1 && regnm(ii,1,kk)==1 && regArcnm(ii,2,kkp)==1 && regnm(ii,2,kkp)==1 && regA(ii,1,kk)>regA(ii,1,kkp) )
            cond2b = 1;
            break;
        end
    end
    %     cond2 = ( cond2a || and(cond2b,condlgmopkspt(ii)) );
    [~,ind1] = max(regA(ii,1,1:2)); %find the larger region
    [~,ind2] = min(regA(ii,1,1:2)); %find the smaller region
    cond2a = ( or(regnm(ii,1,ind1)==1,regnm(ii,2,ind1)==1) && or(regArcnm(ii,1,ind1)==1,regArcnm(ii,2,ind1)==1) ); %iii condition page 182; if there is one polarity region at either the 1st or 2nd thrshold AND there is one arc at either the 1st or 2nd thrshold
    cond2b1 = ( regnm(ii,1,ind2)==1 && regArcnm(ii,1,ind2)==1 ); %Condition ii part 1 from paper
    
    % % changed by Kaylee and Ahmad
    %cond2b2 = ( regnm(ii,1,ind2)==0 && regArcnm(ii,1,ind2)==0 && regY(ii,1,ind1)<=0 ); %Condition ii part 2 from paper
    cond2 = (cond2a && cond2b1);
    %cond2b2 = ( regnm(ii,1,ind2)==0 && regArcnm(ii,1,ind2)==0 && regY(ii,1,ind1)<=0 ); %Condition ii part 2 from paper
    %cond2 = (cond2a && or(cond2b1,cond2b2));
    cond3a = ( regcnt2org(ii,1)<0.3 ); %Condition 3: d <0.3 page 6 supp aEREMCOR
    cond3b = 0;
    if( cond3a==0 && cond2b==1 )
        reg12cnt2org(ii) = lne2pt([regX(ii,1,kk) regY(ii,1,kk) 0],[regX(ii,2,kkp) regY(ii,2,kkp) 0],[0 0 0]); %find distance between these points d'
        if( lt(reg12cnt2org(ii),0.3) ) % d'<0.3
            cond3b = 1;
        end
    end
    cond3 = ( cond3a || and(cond3b,condlgmopkspt(ii)) ); % save locs where either cond3a or (cond3b and the bcg identified ICs from spectrum analysis) are true
    cond4a = ( all(gt(regA(ii,1,1:2),0.1)) && all(gt(regArcA(ii,1,1:2),0.25)) ); %Condition 4: sets the minimum area of the polarity arc region parameters from table S4
    cond4b = ( max(regA(ii,1,1:2))>0.1 && min(regA(ii,1,1:2))>0.05 && max(regArcA(ii,1,1:2))>0.25 && min(regArcA(ii,1,1:2))>0.08 ); %Condition 4: sets the minimum area of the secondary polarity region parameters from table S4
    cond4c = 0;
    if( cond4a==0 && cond4b==0 && cond2b==1 ) % if all of condition 4 is met...
        if( regA(ii,1,kk)>0.1 && regA(ii,2,kkp)>0.1 && regArcA(ii,1,kk)>0.25 && regArcA(ii,2,kkp)>0.12 ) %and if the area of the regions are above a certain percentage
            cond4c = 1; % yes to Condition 4
        end
    end
    cond4 = ( regA(ii,2,ind1)>0.25 && regArcA(ii,2,ind1)>.025 ); %Condition 4 from Cardiac Paper is where these areas are above a certain percentage
    if( cond1 && cond2 && cond3 && cond4 ) %if all conditions are met
        if((cond3a||cond3b) && (cond4a||cond4b||cond4c)) %if one of the 3rd and 4th conditions are met
            condrmmot(ii) = 1; % note that location with a 1
        end
    end
    if( cond1 && cond2 && cond4) %&& ~cond5) %if condition 1, condition 2, and condition 4 are met
        cbtpind(ii) = 1; % cbtpind tells us the IC numbers that are classified as BCG
    end
    
end

% removal of blink and saccade components
condbk = zeros(1,icnum);
condsc = zeros(1,icnum);
condbk(intersect(find(condbkspt==1),find(condbktp==1))) = 1;
mofrqind(condbk==1) = 0;
motpind(condbk==1) = 0;
condlgmopkspt(condbk==1) = 0;
cbtpind(condbk==1) = 0;
condsc(intersect(find(condscspt==1),find(condsctp==1))) = 1;
mofrqind(condsc==1) = 0;
motpind(condsc==1) = 0;
condlgmopkspt(condsc==1) = 0;
cbtpind(condsc==1) = 0;

% signal contribution analysis
%rmpos
rmsignifind = zeros(1,icnum);
cbsignifind = zeros(1,icnum);
signifcond = zeros(icnum,nbch);
spkdurintv = 0.04/mriperd*szsm; %number of data point in 40ms duration
rmpos = zeros(szic(1),szic(2));
rmneg = zeros(szic(1),szic(2));
cbpos = zeros(szic(1),szic(2));
cbneg = zeros(szic(1),szic(2));
for ii=1:icnum
    rmpos(ii,ic(ii,1:szsm)>icave(ii)+4*icstd(ii)) = 1;
    rmneg(ii,ic(ii,1:szsm)<icave(ii)-4*icstd(ii)) = 1;
    cbpos(ii,ic(ii,1:szsm)>icave(ii)+0.1*icstd(ii)) = 1;
    cbneg(ii,ic(ii,1:szsm)<icave(ii)-0.1*icstd(ii)) = 1;
end
eachperd = floor(10*srte);
perdnm = floor(mriperd*srte/eachperd);
maxsig = zeros(nbch,1);
meansig = zeros(nbch,1);
stdsig = zeros(nbch,1);
minstdsig = zeros(nbch,1);
rmsig = zeros(nbch,2);
nrmperd = zeros(nbch,szsm);
for jj=1:nbch
    maxsig(jj) = max(abs(mixsig(jj,1:szsm)));
    meansig(jj) = mean(mixsig(jj,1:szsm));
    stdsig(jj) = std(mixsig(jj,1:szsm));
    minstdsig(jj) = std(mixsig(jj,1:eachperd));
    for kk=2:perdnm
        minstdsig(jj) = min(minstdsig(jj),std(mixsig(jj,(kk-1)*eachperd+1:kk*eachperd)));
    end
    rmsig(jj,1) = meansig(jj) + 4 * stdsig(jj);
    rmsig(jj,2) = meansig(jj) - 4* stdsig(jj);
    nrmperd(jj,mixsig(jj,1:szsm)<=meansig(jj)+4*minstdsig(jj) & mixsig(jj,1:szsm)>=meansig(jj)-4*minstdsig(jj) ) = 1;
end
flgSingleChArt=zeros(size(ic,1),1);
for ii = 1:icnum
    [possectnm,possectrge] = cont1(rmpos(ii,1:szsm));
    [negsectnm,negsectrge] = cont1(rmneg(ii,1:szsm));
    selic = zeros(size(ic));
    selic(ii,1:szsm) = ic(ii,1:szsm);
    corrsig = mixsig - A*selic;
    
    % added by ahmad and kaylee:
    subsig=A*selic;
    [spt1,frq1] = pwelchcal(subsig,srte,512,256); %spt is the power spectral density curve, frq is the freq. of the PSD curve
    spt1_sum = sum(spt1,2);
    [spt1sum spt1sumind] = sort(spt1_sum,'descend');
    spt1maxes = spt1sum(1:3);

    if (spt1maxes(1)>5*spt1maxes(2)||spt1maxes(1)>10*spt1maxes(3)) && ickurt(ii)>4  && alpha_nm<max(max(delta_nm,theta_nm),beta_nm)
        
        flgSingleChArt(ii)=1;
    end

    
    RapidMotionCondition = ( ickurt(ii)>5 && (sum(rmpos(ii,1:szsm))>0 || sum(rmneg(ii,1:szsm))>0) );
    if( RapidMotionCondition )
        cond1 = zeros(nbch,possectnm);
        for kk=1:possectnm
            if(possectrge(kk,2)-possectrge(kk,1)>spkdurintv) % check if the duration is more than 40 ms
                Diff = max(abs(mixsig(:,possectrge(kk,1):possectrge(kk,2))-corrsig(:,possectrge(kk,1):possectrge(kk,2))),[],2); %delta_jk; The signal reduction at the iith IC
                RelDiff = Diff ./ max(abs(mixsig(:,possectrge(kk,1):possectrge(kk,2))),[],2);%delta_jk/Vjk,max
                AbsDiff = Diff ./ max(maxsig(:),100*ones(size(maxsig))); %????
                Ave1 = abs(mean(mixsig(:,possectrge(kk,1):possectrge(kk,2)),2));
                Ave2 = abs(mean(corrsig(:,possectrge(kk,1):possectrge(kk,2)),2));
                SignalSpike = or(gt(max(mixsig(:,possectrge(kk,1):possectrge(kk,2)),[],2),rmsig(:,1)),lt(min(mixsig(:,possectrge(kk,1):possectrge(kk,2)),[],2),rmsig(:,2)));
                Significance = ( ( lt(Ave2,0.9*Ave1) & gt(RelDiff,diffr(1)*ones(size(RelDiff))) ) & ( gt(AbsDiff,diffr(2)*ones(size(AbsDiff))) | gt(Diff,200*ones(size(Diff))) ) );
                cond1(:,kk) = and( Significance(:), SignalSpike(:) );
                for jj=1:nbch
                    if(cond1(jj,kk)==1)
                    end
                end
            end
        end
        cond1(cond1(:, sum(cond1(:,:),1) < 2) == 1) = 0;
        signifcond(ii,:) = ( sum(cond1,2)>0 );
        
        cond1 = zeros(nbch,negsectnm);
        for kk=1:negsectnm
            if(negsectrge(kk,2)-negsectrge(kk,1)>spkdurintv)
                Diff = max(abs(mixsig(:,negsectrge(kk,1):negsectrge(kk,2))-corrsig(:,negsectrge(kk,1):negsectrge(kk,2))),[],2);
                RelDiff = Diff ./ max(abs(mixsig(:,negsectrge(kk,1):negsectrge(kk,2))),[],2);
                AbsDiff = Diff ./ max(maxsig(:),100*ones(size(maxsig)));
                Ave1 = abs(mean(mixsig(:,negsectrge(kk,1):negsectrge(kk,2)),2));
                Ave2 = abs(mean(corrsig(:,negsectrge(kk,1):negsectrge(kk,2)),2));
                SignalSpike = or(gt(max(mixsig(:,negsectrge(kk,1):negsectrge(kk,2)),[],2),rmsig(:,1)),lt(min(mixsig(:,negsectrge(kk,1):negsectrge(kk,2)),[],2),rmsig(:,2)));
                Significance = ( ( lt(Ave2,0.9*Ave1) & gt(RelDiff,diffr(1)*ones(size(RelDiff))) ) & ( gt(AbsDiff,diffr(2)*ones(size(AbsDiff))) | gt(Diff,200*ones(size(Diff))) ) );
                cond1(:,kk) = and( Significance(:), SignalSpike(:) );
                for jj=1:nbch
                    if(cond1(jj,kk)==1)
                    end
                end
            end
        end
        cond1(cond1(:, sum(cond1(:,:),1) < 2) == 1) = 0;
        signifcond(ii,:) = or( signifcond(ii,:), transpose(sum(cond1,2)>0) );
    end
    indrm = ( signifcond(ii,ignsignif==0)==1 );
    if( condrmmot(ii) && sum(indrm)>=6 )
        rmsignifind(ii) = 1;
    end
end

% bcg ic
cond3 = zeros(icnum,nbch);
for ii=1:icnum
    selic = zeros(size(ic));
    selic(ii,1:szsm) = ic(ii,1:szsm);
    corrsig = mixsig - A*selic; % Corrected Signal after removing IC ii
    if( sum(cbpos(ii,1:szsm))>0 && sum(cbneg(ii,1:szsm))>0 )
        for jj=1:nbch
            posave1 = mean(abs( mixsig(jj,cbpos(ii,1:szsm)==1 & nrmperd(jj,1:szsm)==1) ));
            posave2 = mean(abs( corrsig(jj,cbpos(ii,1:szsm)==1 & nrmperd(jj,1:szsm)==1) ));
            negave1 = mean(abs( mixsig(jj,cbneg(ii,1:szsm)==1 & nrmperd(jj,1:szsm)==1) ));
            negave2 = mean(abs( corrsig(jj,cbneg(ii,1:szsm)==1 & nrmperd(jj,1:szsm)==1) ));
            cond2 = ( min(posave2/posave1,negave2/negave1)<0.95 && 0.5*(posave2/posave1+negave2/negave1)<0.97 );
            if(cond2)
                signifcond(ii,jj) = signifcond(ii,jj) + 2;
            end
            cond3(ii,jj) = ( 0.5*(posave2/posave1+negave2/negave1)<0.8 );
        end
        indcb = ( signifcond(ii,:)==2 | signifcond(ii,:)==3 );
        if( condlgmopkspt(ii) && sum(indcb)>=1 )
            cbsignifind(ii) = 1;
        end
    end
    signifcond(ii,find(cond3(ii,:)==1)) = signifcond(ii,find(cond3(ii,:)==1)) + 4;
    indrsp = ( signifcond(ii,ignsignif==0)>=4 );
end

% selection combination
cbicind = ( condlgmopkspt==1 & cbtpind==1 & cbsignifind==1 );
%final condition for spectrum & the topographic map
condkham = (cond5 && (alpha_nm > max(max(delta_nm,theta_nm),beta_nm))) || (cond5 && ((alphaband(1) < MaxPkLoc)&& (alphaband(2) > MaxPkLoc)));
for j=1:size(cbicind,2)
    if condkham == 1
        cbicind(j) = 0;
    end
end
topomap = cbtpind;
spectra_BCG_ICs = condlgmopkspt;
blink_IC = condbk;
saccade_IC = condsc;
smolregion = toosmall;
tpblink = condbktp;
tpsac = condsctp;
singchan=flgSingleChArt;
muscleloc=gammaflg;
end
