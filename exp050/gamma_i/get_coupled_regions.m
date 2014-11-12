function regions=get_coupled_regions(A)

[m,p]=size(A);
L=zeros(p,1); % label vector
wet=true(p,1); % all 
iregion=1; % region label

while 1
    ii=find(wet==1,1,'first');
    if isempty(ii); break; end
    idx=ii; % linear indices of the pixels that were just found
    wet(idx)=0; % set the newly found pixels to dry
    L(idx) = iregion; % label first point
    n=[];
    while ~isempty(idx); % find neighbours

        cols=A(:,idx)~=0;
        rows=any(cols,2);
        hoit=A(rows,:)~=0;

        n=find(any(hoit,1));
        
        n = unique(n(:)); % remove duplicates
        idx = n(wet(n)); % keep only the neighbors that haven't been labelled
        if isempty(idx); break; end
        L(idx) = iregion;
        wet(idx)=0; % set the pixels that were just labeled to 0

    end
    iregion=iregion+1;
end

for ii=1:iregion-1;
    regions{ii}=find(L==ii);
end

