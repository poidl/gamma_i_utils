function cregs=find_coupled(A)
    %error('this doesn''t work I think?')
    [m,n]=size(A);
    cregs=nan*ones(n,1);

    done_r=[];
    list_r=1; % row list
    ireg=1;
    %keyboard
    while any(isnan(cregs))
        while ~isempty(list_r);
            r=list_r(1);
            cols=find(A(r,:)~=0);
            cregs(cols)=ireg;
            sum(~isnan(cregs))
            for c=cols
                rows=find(A(:,c)~=0);
                list_r=[list_r;rows];
            end   
            list_r=unique(list_r);            
            done_r=[done_r,r];
            list_r=setdiff(list_r,done_r);
        end
        %keyboard
        A=A(~done_r,:);
        done_r=[];
        list_r=1; % row list
        ireg=ireg+1;
    end 
end
