% suspected motion segment
function bd2ind = badsegmo(EEG)
sdeind = ones(1,EEG.nbchan);
sdeind([17:19,21:24,32]) = 0;
Index = EEG.bad(find(sdeind==1),:);
Index(Index(:, sum(Index(:,:),1) < 2)==1) = 0;
bd2ind = ( sum(Index,1)>=6 );
end


