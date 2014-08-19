% plotting script
clear all;
close all;

%'gins3d','pns3d','p_s','gi_bb','p_bb','ilat','ilon'
addpath(genpath('../../../../omega/ansu_utils/external_scripts/'))
load('../data/gins3d.mat')
load('../data/gamma_initial.mat')

[ns,ny,nx]=size(pns3d);
dp3d=nan*pns3d;
for kk=1:ns
   dp3d(kk,:,:)= p_s(kk,:,:)-pns3d(kk,:,:);
end 

va=dp3d;
%va=va(ns-7:ns,:,:);
lat=squeeze(lat(1,:,1));
lon=squeeze(lon(1,1,:));

%va=va(1:6,:,:);
%va=va(ns-6:ns,:,:);
vv=va; % variable to plot
%vv=diff(vv,1);
nit=size(vv,1);

nfig=floor(nit/6)+(mod(nit,6)>0); % number of figures (pages)
ncols=2; % number of columns
nrows=3; % number of rows
nsp=ncols*nrows; % max. number of subplots per figure
%nsp=1

spwidth=0.4; % subplot width
spheight=0.2; % subplot height
wscol= 0.08; % white space between columns
wsrow=0.1; % white space between rows
leftmarg=(1-(ncols*spwidth+(ncols-1)*wscol))*0.5; % left and right margin
topmarg=(1-(nrows*spheight+(nrows-1)*wsrow))*0.5; % top and bottom margin

cmax=log10(max(abs(vv(:))));
%cmax=-5
cmin=-13;
cbh=nan*ones(nsp); % colorbar handles
fac=nan*ones(nsp);
cmax2=max(abs(vv(:)));
%cmax2=1e-4



colorscale={'lin'};
txt={'a)','b)','c)','d)','e)','f)'};
for icolorscale=1:1
    iit=1; % index of iteration
    
    for ifig=1:nfig
        sz=2.0*[21 29.7];
        figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)],'Visible','off') 
        set(gcf,'DefaultAxesFontSize', 15)
        set(gcf,'DefaultTextFontSize',15)

        isp=1; % index of subplot
        ip=1; % index of plot
        
        cmp=colormap(hot(128));
        cmp2=cmp(find(cmp(:,1)==1,1,'first')-15:end,:);
        
        while (isp<=nsp) && ( ip <=nit) && (iit<=size(vv,1))
            ip=((ifig-1)*nsp+isp);
            irow=ceil(isp/ncols); % row index
            icol=isp-(irow-1)*ncols; % column index
            left=leftmarg+(icol-1)*(spwidth+wscol); % current subplot position
            bottom=1-topmarg-irow*spheight-(irow-1)*wsrow; % current subplot position

            subplot('position',[left,bottom,spwidth,spheight])
            vp=squeeze(vv(iit,:,:));

            if colorscale{icolorscale}=='log'
                tag='log';
                tmp=vp;
                tmp(vp<=0)=nan;
                pos=log10(   tmp );
                h=imagesc(lon,lat,pos);
                set(h,'alphadata',~isnan(pos)) % white nans
                set(gca,'YDir','normal')
                colormap(flipud(colormap(cmp2))) ;
                caxis([cmin cmax])
                cb=colorbar('location','southoutside','position',[left,bottom-0.35*wsrow,spwidth,0.1*wsrow],'XTickLabel',[]);
                title(['Iteration ',num2str(ip)])
                cbfreeze(cb)
                freezeColors
                hold on

                vp(vp>=0)=nan;
                vp=log10(  -vp );
                h=imagesc(lon,lat,vp);
                colormap(rot90(colormap(cmp2),2)) ;
                caxis([cmin cmax])
                cb=colorbar('location','southoutside','position',[left,bottom-0.5*wsrow,spwidth,0.1*wsrow]);
                xlabel(cb,'Red: log10($\Phi''>0$)  $\phantom{xxxxxx}$ Blue: log10( -1$\cdot(\Phi''<0$))','interpreter','latex','fontsize',18)
                
            elseif colorscale{icolorscale}=='lin'
                tag='lin';

                fac(isp)=max(abs(vp(:))); vp=vp/fac(isp); 
                %vp=vp/cmax2;

                h=imagesc(lon,lat,vp);
                colormap([fliplr(cmp2);flipud(cmp2)]) ;
                caxis([-1 1])
                cbh(isp)=colorbar('Location','SouthOutside','position',[left,bottom-0.4*wsrow,spwidth,0.1*wsrow]);
                if (isp==6) || (iit==size(vv,1)) % 
                    for ii=1:nsp
                        set(cbh(ii),'XTick',-1:0.25:1)
                        %set(cbh(ii),'XTickLabel',num2str(fac(ii)*str2num(get(cbh(ii),'XTickLabel')),'%2.1e'));
                        set(cbh(ii),'XTickLabel',num2str(fac(ii)*str2num(get(cbh(ii),'XTickLabel')),'%0.1e'));
                        %set(cbh(ii),'XTickLabel',num2str(cmax2*str2num(get(cbh(ii),'XTickLabel')),'%2.1e'));
                        %xlabel(cbh(ii),'$\Phi''$','interpreter','latex','fontsize',18)
                        %xlabel(cbh(ii),['p_{ns}(it=',num2str(ip-nsp+ii-1),')-p_{ns}(final)  [db]'],'fontsize',18)
                        %xlabel(cbh(ii),['depth at backbone: ',num2str(pns3d(ip-nsp+ii,ilat,ilon)),'  [dbar]'],'fontsize',18)
                    end
                end
                hold on
                % plot backbone
                plot(lon(ilon),lat(ilat),'k*','markersize',20)
            end

            set(h,'alphadata',~isnan(vp)) % white nans
            set(gca,'YDir','normal')
             load ('coast_data.mat');
             plot(coast_data_long,coast_data_lat,'k-','LineWidth',1);
            title(['pressure difference to neutr. surf. (surface ',num2str(ip),')'])
            xlabel(['depth at backbone: ',num2str(pns3d(ip,ilat,ilon)),'  [dbar]'],'fontsize',18)
            text(-20,95,txt(isp),'fontsize',18')

            iit=iit+1; isp=isp+1;
        end

        print('-dpng','-r200',['../figures/dp3d/',num2str(ifig,'%02i')])
    end
end
