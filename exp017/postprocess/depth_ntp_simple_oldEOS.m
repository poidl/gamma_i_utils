function [sns,ctns,pns] = depth_ntp_simple(s0,ct0,p0,s,ct,p,drho)

% warning('no check of input dimensions')

% search in one direction only, based on the initial value of F (at bottle
% depth). Note that this function will generally not find all zero
% crossings (possibly including stable ones), and may find unstable zero crossings.

%user_input; % get delta
delta=1e-11;

nxy=size(s,2);

inds=1:nxy;

pns = nan(1,nxy);
sns = nan(1,nxy);
ctns = nan(1,nxy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% discard land; inds contains the indices of remaining points
[s0,ct0,p0,s,ct,p,drho,inds]=discard_land(s0,ct0,p0,s,ct,p,drho,inds);

nz=size(s,1);
nxy=size(s,2);
inds_wet=1:nxy;

p0_stacked=repmat(p0(:)',[nz 1]);

k=sum(p0_stacked>=p,1); % adjacent bottle looking up. values 1 to nz.
k3d=k+nz*(0:nxy-1); % 3-d
kinc=nan*ones(1,nxy); % increment of k

p_=p(k3d);
s_=s(k3d);
ct_=ct(k3d);

go_shallow_old=false(1,nxy);
go_deep_old=false(1,nxy);

bisect=false(1,nxy);
done=false(1,nxy);
cnt=0;
while any(~done)
    cnt=cnt+1;
    pmid=0.5*(p0+p_); 
    bottle = rho_from_theta(s0,ct0,pmid);
    F=rho_from_theta(s_,ct_,pmid)-bottle+drho;

    go_shallow = F >=  delta; % go shallower to find first bottle pair
    go_deep    = F <= -delta; % go deeper to find first bottle pair
    final= abs(F)<delta; % hit zero
    
    sns(inds(inds_wet(final)))=s_(final);
    ctns(inds(inds_wet(final)))=ct_(final);
    pns(inds(inds_wet(final)))=p_(final);

    cfb = (go_shallow_old & go_deep); % crossed from below
    cfa = (go_deep_old & go_shallow); % crossed from above
    crossed= cfb|cfa;
    start_bis = (crossed & ~bisect); % start bisection here
    bisect = (bisect | start_bis); % bisect here
    
    search_initial = (go_deep|go_shallow) & ~bisect & ~final;

    kinc(inds_wet(go_deep & search_initial))= 1;
    kinc(inds_wet(go_shallow & search_initial))=-1;
    
    
    k=k+kinc;
    k3d=k3d+kinc;
    k_=k(inds_wet(search_initial));%+kinc(inds_wet(search_initial));
    k3d_=k3d(inds_wet(search_initial));%+kinc(inds_wet(search_initial)); 

    out=(k_<1)|(k_>nz); % above or below domain
    out2=false(1,length(final));
    out2(search_initial)= out;

    done=isnan(F)|final|out2; 
    
    k3d_=k3d_(~out);    
    search_initial=search_initial(~done);
    bisect=bisect(~done);
    crossed=crossed(~done);

    
    if ~all((search_initial|bisect))
        error('something is wrong')
    end  
    
    
    if cnt>1
        s1=s1(~done);
        ct1=ct1(~done);
        p1=p1(~done);
        s2=s2(~done);
        ct2=ct2(~done);
        p2=p2(~done);
        s2(crossed)=s1(crossed);
        ct2(crossed)=ct1(crossed);
        p2(crossed)=p1(crossed);
    else
        s2=s_(~done);
        ct2=ct_(~done);
        p2=p_(~done);
    end
    s1=s_(~done);
    ct1=ct_(~done);
    p1=p_(~done);   

%%%%%%%%%%%%%%%%%%
% data for next evaluation of F
    s0=s0(~done);
    ct0=ct0(~done);
    p0=p0(~done);

    drho=drho(~done);
    s_=nan*ones(1,sum(~done));
    ct_=nan*ones(1,sum(~done));
    p_=nan*ones(1,sum(~done));

    s_(search_initial)=s(k3d_);
    ct_(search_initial)=ct(k3d_);
    p_(search_initial)=p(k3d_); 
    
    if any(bisect)
        s_(bisect) =0.5*(s1(bisect) +s2(bisect));
        ct_(bisect)=0.5*(ct1(bisect)+ct2(bisect));
        p_(bisect) =0.5*(p1(bisect) +p2(bisect));
    end
%%%%%%%%%%%%%%%%%%    
 
    go_shallow_old=go_shallow(~done);
    go_deep_old=go_deep(~done);
    
    inds_wet=inds_wet(~done);   
end

end




function [s0,ct0,p0,s,ct,p,drho,inds]=discard_land(s0,ct0,p0,s,ct,p,drho,inds)
    iwet=~isnan(s0);

    s0=s0(iwet);
    ct0=ct0(iwet);
    p0=p0(iwet);
    drho=drho(iwet);
    s=s(:,iwet);
    ct=ct(:,iwet);
    p=p(:,iwet);
    inds=inds(iwet);
end
