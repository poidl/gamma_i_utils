rpot=gsw_rho(s,ct,0*p);

lat=squeeze(lats(1,:,1));

pp=squeeze(p(:,1,1));


va=rpot;
vals=plt_get_values(va(:),5,10);
vals_label=num2str(vals.','%.2f');
vals=str2num(vals_label);

va=squeeze(va(:,:,1));
%h=imagesc(lats(1,:,1),pp(:,1,1),va);
h=imagesc(lat,pp,va);
set(h,'alphadata',~isnan(va))
hold on 
[c,h]=contour(lat,pp,va,vals,'color','k');
clabel(c,h,'labelspacing',1e10)

colorbar()
title('\sigma_0')
print('-dpng','-r200',['figures/rpot.png'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
grf=gamma_rf(s,ct);

va=grf;
vals=plt_get_values(va(:),5,10);
vals_label=num2str(vals.','%.2f');
vals=str2num(vals_label);

va=squeeze(va(:,:,1));
%h=imagesc(lats(1,:,1),pp(:,1,1),va);
h=imagesc(lat,pp,va);
set(h,'alphadata',~isnan(va))
hold on 
[c,h]=contour(lat,pp,va,vals,'color','k');
clabel(c,h,'labelspacing',1e10)

colorbar()
title('\gamma_{rf}')

print('-dpng','-r200',['figures/grf.png'])
