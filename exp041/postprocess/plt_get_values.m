function values=plt_get_values(va,nbins,nsurf)
    % get vector of values of va whose spacing decreases with
    % stratification (decreses with number of observations in a density range)

    va_min=min(va(:));
    va_max=max(va(:))+eps(1100); % add eps(1100) such that histc() returns zero for last bin.
    bins=linspace(va_min,va_max,nbins);
    hi=histc(va(:),bins);
    hi=hi(1:end-1); % remove last entry
    %nsurf=30; % approx. number of total surfaces
    N=int16(hi*nsurf/sum(hi)); % number of surfaces for bin
    N=double(N);
    dg=diff(bins)./N'; % increment of gamma
    dg(~isfinite(dg))=nan;
    
    nsurf=sum(N);
    values=nan*(ones(1,nsurf));
    
    jj=1;
    for ii=1:length(dg);
        if ~isnan(dg(ii))
            values(jj: jj+N(ii)-1)=bins(ii)+dg(ii)*(0.5 : 1 : N(ii)-0.5);
            jj=jj+N(ii);
        end
    end 
end
