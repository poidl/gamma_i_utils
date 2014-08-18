close all
clear all

sz=1.5*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 

[s,ct,p]=gammanc_to_sctp();
%load('/home/nfs/z3439823/mymatlab/paul/paul_data/WGHC/wghc_2deg_2010.mat')
%load('/home/nfs/z3439823/mymatlab/paul/paul_data/WGHC/wghc_2degree.mat')

% for pdf file names
data_set_string='_gammanc'
% set limit for close-up plot
limit=1e-6;

[nz,ny,nx]=size(s);

[n2,pmid]=gsw_Nsquared(s(:,:),ct(:,:),p(:,:));
n2=reshape(n2,[nz-1,ny,nx]);
pmid=reshape(pmid,[nz-1,ny,nx]);


ip=n2>=0;
in=n2<0;
n2p=n2(ip);
n2n=n2(in);

n2_tiny=abs(n2(:))<limit;
n2t=n2(n2_tiny);


xl1=min(n2(:));
xl2=max(n2(:));
x=linspace(xl1,xl2,100);
dx=diff(x); dx=dx(1);
x=x-min(x(x>0));
x=[x,x(end)+dx];

xp=x(x>=0);
xn=x(x<=0); 
[np,vals]=histc(n2p(:),xp); % histc: cut last element and shift x-axis for plot
np=np(1:end-1); xp=xp(1:end-1)+0.5*diff(xp);
[nn,vals]=histc(n2n(:),xn);
nn=nn(1:end-1); xn=xn(1:end-1)+0.5*diff(xn);

n=[nn;np];
x=[xn,xp];


bv=0.5;
hb=bar(x, n,  'barwidth', 1,'basevalue', bv)

% color the negative values in red
ch = get(hb,'Children');
fvd = get(ch,'Faces');
vvd = get(ch,'Vertices'); % one indexed color per vertex: http://www.mathworks.com.au/help/matlab/ref/patch_props.html#Faces
fvcd = get(ch,'FaceVertexCData');
fvcd(fvd(1:length(xn),:))=1.8; % the color axis should range from 1 to 2 (check with caxis)
set(ch,'FaceVertexCData',fvcd);

set(gca,'YScale','log')
xlim([xl1-0.1*xl2,1.1*xl2]);
ylim([log(bv),1.1*max(np)])
xlabel('N^2')
ylabel('frequency')
% subplot(2,1,2)
% bar(xn, -nn, 'barwidth', 1, 'basevalue', bv);
% %set(gca,'YScale','log')
% xlim([xl1-0.1*xl2,1.1*xl2]);
% ylim([log(bv),1.1*max(np)])


print('-dpdf','-r200',['figures/n2_hist',data_set_string,'.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sz=1.5*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)])

%xl1=min(min(n2(:),-limit));
xl1=max(min(n2(:),-limit));
xl2=min(max(n2(:),limit));
x=linspace(xl1,xl2,100);
dx=diff(x); dx=dx(1);
x=x-min(x(x>0));
x=[x,x(end)+dx];

xp=x(x>=0);
xn=x(x<=0); 
[np,vals]=histc(n2p(:),xp); % histc: cut last element and shift x-axis for plot
np=np(1:end-1); xp=xp(1:end-1)+0.5*diff(xp);
[nn,vals]=histc(n2n(:),xn);
nn=nn(1:end-1); xn=xn(1:end-1)+0.5*diff(xn);

n=[nn;np];
x=[xn,xp];


bv=0.5;
hb=bar(x, n,  'barwidth', 1,'basevalue', bv)

% color the negative values in red
ch = get(hb,'Children');
fvd = get(ch,'Faces');
vvd = get(ch,'Vertices'); % one indexed color per vertex: http://www.mathworks.com.au/help/matlab/ref/patch_props.html#Faces
fvcd = get(ch,'FaceVertexCData');
fvcd(fvd(1:length(xn),:))=1.8; % the color axis should range from 1 to 2 (check with caxis)
set(ch,'FaceVertexCData',fvcd);

set(gca,'YScale','log')
xlim([xl1-0.1*xl2,1.1*xl2]);
ylim([log(bv),1.1*max(np)])

xlabel('N^2')
ylabel('frequency')


xp1=x(n>0 & x'>0);
xn1=x(n>0 & x'<0);

yl=get(gca,'ylim');
hl1=line([xp1(1),xp1(1)]-0.5*dx,[bv,yl(2)],'color','b')
hl2=line([xn1(end),xn1(end)]+0.5*dx,[bv,yl(2)],'color','r')
str1=['x=',num2str(xn1(1)+0.5*dx)];
str2=['x=',num2str(xp1(1)-0.5*dx)];
str1=['min(N^2>=0)=',num2str(min(n2p))];
str2=['max(N^2<0)=',num2str(max(n2n))];
legend([hl1,hl2],str1,str2,'location','northwest')

axis tight
print('-dpdf','-r200',['figures/n2_hist_tiny',data_set_string,'.pdf'])








