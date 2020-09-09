
function [pkind,sel,selflg,Bp2] = selpk2(Bp2,itvl,coritvl,coritvlpkloc,thr,cond)
pkind = find(((Bp2.pkloc>=itvl(1) & Bp2.pkloc<itvl(2)))==1);
if(numel(pkind)==0)
    sel = [];
    return;
end
cond(3) = (cond(3) && numel(coritvlpkloc)>0);
pkloc = Bp2.pkloc(pkind);
selflg = zeros(numel(pkloc),numel(cond));
if( cond(1)>0 )
    [~,Ind] = max(max([Bp2.lfslp(pkind);Bp2.rgslp(pkind)],[],1));
    selflg(Ind,1) = 1;
end
if( cond(2)>0 )
    [~,Ind] = max(Bp2.pkamp(pkind));
    selflg(Ind,2) = 1;
end
if( cond(3)>0 )
    maxcor = zeros(size(coritvlpkloc));
    for jj=1:size(coritvl,1)
        [Corr,lag] = xcorr(Bp2.data(coritvl(jj,1):coritvl(jj,2)),Bp2.data(itvl(1):itvl(2)));
        [~,ind] = max(Corr);
        maxcor(jj) = itvl(1)+coritvlpkloc(jj)-coritvl(jj,1)-lag(ind);
    end
    selsc = abs(pkloc-mean(maxcor));
    [~,Ind] = min(selsc);
    selflg(Ind,3) = 1;
end
if( cond(4)>0 )
    for ii=1:numel(pkloc)
        selflg(ii,4) = ( Bp2.pkamp(pkind(ii))>=thr(1) || max(Bp2.lfslp(pkind(ii)),Bp2.rgslp(pkind(ii)))>=thr(2) );
    end
else
    selflg(:,4) = 1;
end

sel = ones(numel(pkloc),1);
for ii=1:numel(cond)
    if(cond(ii)>0)
        sel = ( sel & selflg(:,ii) );
    end
end
if(numel(find(sel==1))==0)
    if(cond(3)==1)
        sel = ( selflg(:,1) & selflg(:,2) & selflg(:,3) );
        if(numel(find(sel==1))==0)
            sel = ( selflg(:,4) & selflg(:,3) );
            if(numel(find(sel==1))==0)
                ind1 = [find(selflg(:,2)&selflg(:,3)==1),find(selflg(:,1)&selflg(:,3)==1)];
                if(numel(ind1)==0)
                    ind1 = [find(selflg(:,2)==1),find(selflg(:,1)==1)];
                end
                [~,ind] = min(selsc(ind1));
                sel(ind1(ind)) = 1;
            end
        end
    elseif(cond(3)==0)
        sel = ( selflg(:,1) & selflg(:,2) & selflg(:,4) );
        if(numel(find(sel==1))==0)
            ind1 = find(selflg(:,2)&selflg(:,4)==1);
            if(numel(ind1)==0)
                ind1 = find(selflg(:,1)&selflg(:,4)==1);
                if(numel(ind1)==0)
                    ind1 = find(selflg(:,2)==1);
                end
            end
            sel(ind1) = 1;
        end
    end
end
end
