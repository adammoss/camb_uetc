function plot_uetc_recon(filename,type,norm,titlename)

lab_fontsize =12; axes_fontsize =12;
nn = 256;
clf;
%filename = 'v05_n256_m32_w025'
%type = 'vv'

ktau = importdata(strcat('data/',filename,'_uetc_ktau.dat'));
n = numel(ktau);

emt = importdata(strcat('data/',filename,'_uetc_',type,'.dat'));
emt_recon = importdata(strcat('data/',filename,'_uetc_',type,'_recon.dat'));
emt_diff = emt - emt_recon;

subplot(2,2,1);

log_ktau = log10(ktau); 
xlin = linspace(min(log_ktau),max(log_ktau),nn);
[X,Y] = meshgrid(xlin,xlin);
Z = griddata(log_ktau,log_ktau,emt,X,Y,'cubic');
[ch,ch]=contourf(X,Y,Z,400);
set(ch,'edgecolor','none');
xlabel('log_{10} (k \tau_1)','FontSize',lab_fontsize);
ylabel('log_{10} (k \tau_2)','FontSize',lab_fontsize);
%zlabel('(\tau_1 \tau_2)^{1/2} <\theta^V (k \tau_1) \theta^V (k \tau_2)> ','FontSize',lab_fontsize);
%axis([-2,2.3,-2,2.3])
axis([1.5,2.6,1.5,2.6])
colorbar;
caxis([-2 2])

subplot(2,2,2);

log_ktau = log10(ktau); 
xlin = linspace(min(log_ktau),max(log_ktau),nn);
[X,Y] = meshgrid(xlin,xlin);
Z = griddata(log_ktau,log_ktau,emt_recon,X,Y,'cubic');
[ch,ch]=contourf(X,Y,Z,200);
set(ch,'edgecolor','none');
xlabel('log_{10} (k \tau_1)','FontSize',lab_fontsize);
ylabel('log_{10} (k \tau_2)','FontSize',lab_fontsize);
%zlabel('(\tau_1 \tau_2)^{1/2} <\theta^V (k \tau_1) \theta^V (k \tau_2)> ','FontSize',lab_fontsize);
%axis([-2,2.3,-2,2.3])
axis([1.5,2.6,1.5,2.6])
colorbar;
caxis([-2 2])

set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperPosition',[ 0 0 6 5]);
outfile = strcat('plots/',filename,'_uetc_',type,'.eps');
print('-depsc2',outfile);
