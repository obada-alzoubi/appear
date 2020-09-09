% Independent component analysis
function [A,W, icaweights,icasphere] = ICA(Input,rmbad)
    tICA = tic;
    fprintf('Independent component analysis ...\n');
    if(strcmp(rmbad,'y'))
        Input.data = Input.data(1:31,Input.badmot==0);
    end
    Input = pop_runica(Input, 'icatype', 'runica', 'extended',1,...
        'interrupt','off');
    icaweights =Input.icaweights;
    icasphere =Input.icasphere;
    W = Input.icaweights*Input.icasphere;
    %cmps = W * Input.data;
    A = inv(W);
    toc(tICA);
end
