
% cardiac cycle determination
function [pkloc, mbpfn, A, W,icaweights,icasphere] = cardperiod(input1,acqfrq,ECG_ind)
    tic;
    fprintf('Determining cardiac period ...\n');
    input1.data = input1.data(setdiff(1:end, ECG_ind),:); %Shouold be changed for more than 31 channels 
    [A,W,icaweights,icasphere ] = ICA(input1,'n');
    cbind = icid(W*input1.data,A,input1.data,input1.srate,input1.times(end));

    mbpfn = icsel(W*input1.data,input1.srate,cbind,acqfrq);
    %BCG_filename=sprintf('%s_BCG.csv',EXP);
    %csvwrite(BCG_filename,mbpfn)
    pkloc = mbp(mbpfn,input1.srate,acqfrq);
    toc;
end
