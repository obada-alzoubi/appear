function EEG = preproc(EEG)
    tr =EEG.APEAR.TR;
    slpertr = EEG.APEAR.slice_per_TR;
    mrkA    = EEG.APEAR.mrkA;
    dwnsmr = EEG.APEAR.Fs;
    bpfrq = EEG.APEAR.filterRange;
    slfrq = slpertr / tr;
    temp = [60, 26, slfrq*(1:1:100)];
    temp = temp(temp<=max(2*bpfrq(2),temp(1)));
    temp = sort(temp);
    bsfrq = zeros(numel(temp),2);
        for ii=1:numel(temp)
        bsfrq(ii,1) = temp(ii)-0.5;
        bsfrq(ii,2) = temp(ii)+0.5;
        end %end for loop
    nbch = size(EEG.data,1);

    % gradient artifact removal (using FMRIB plugin)
    tic;
    fprintf('Correcting MRI artifact ...\n');
    EEG = pop_fmrib_fastr(EEG,bpfrq(2),1,30,'R128',0,0,0,(mrkA),slpertr,0.03,32,'auto');
    EEG.pnts = mrkA(end)+round(tr*EEG.srate)-mrkA(1);
    EEG.xmax = EEG.pnts;
    EEG.data = EEG.data(:,mrkA(1):mrkA(1)+EEG.pnts-1);
    EEG.times = (1:1:EEG.pnts)/EEG.srate;
    toc;

    % downsampling to 250S/s
    tic;
    fprintf('Downsampling the eeg data ...\n');
    redsmsz = EEG.pnts*dwnsmr/EEG.srate;
    sgn1 = zeros(nbch,redsmsz);
        for ii=1:nbch
        sgn1(ii,1:redsmsz) = resample(EEG.data(ii,1:EEG.pnts),dwnsmr,EEG.srate,10);
        end %end for loop
    EEG.data = sgn1;
    EEG.srate = dwnsmr;
    EEG.pnts = redsmsz;
    EEG.xmax = EEG.pnts;
    EEG.times = (1:1:EEG.pnts)/EEG.srate;
    clear ReducedSignal;
    toc;

    % bandpass filtering
    tic;
    fprintf('Band-pass frequency filtering the eeg data ...\n');
    bpfrq(2) = min(bpfrq(2),0.5*EEG.srate);
    EEG.data = eegfilt(EEG.data,EEG.srate,bpfrq(1),0,0,3*fix(EEG.srate/bpfrq(1)),0,'fir1',0);
    EEG.data = eegfilt(EEG.data,EEG.srate,0,bpfrq(2),0,3*fix(EEG.srate/bpfrq(1)),0,'fir1',0);
    toc;

    % bandstop filtering at multiple acquisition frequencies, 26Hz, 60Hz
    tic;
    bsfrq = bsfrq(bsfrq(:,2)<0.5*EEG.srate,1:2);
    bsfrq2 = size(bsfrq);
    fprintf('Band-stop frequency filtering the eeg data ...\n');
    for iibs = 1:bsfrq2(1)
        Wn = bsfrq(iibs,1:2)/(0.5*EEG.srate);
        [nu,de] = butter(3,Wn,'stop');
            for ii=1:nbch
                EEG.data(ii,:) = filtfilt(nu,de,EEG.data(ii,:));
            end %end for loop
    end %end for loop
    toc;
end %end preproc function