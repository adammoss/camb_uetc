function plot_uetc(filename,type,norm,titlename)

lab_fontsize =16; axes_fontsize =16; title_fontsize=16;
clf;
%filename = 'v08_n256_m64_w025'
%type = 'vv'

ktau = importdata(strcat('data/',filename,'_uetc_ktau.dat'));
n = numel(ktau);

emt = importdata(strcat('data/',filename,'_uetc_',type,'.dat'));
emt_max = max(max(emt))
if norm==1
  emt=emt/emt_max;
end

log_ktau = log10(ktau); 
xlin = linspace(min(log_ktau),max(log_ktau),1024);
[X,Y] = meshgrid(xlin,xlin);
Z = griddata(log_ktau,log_ktau,emt,X,Y,'cubic');
[ch,ch]=contourf(X,Y,Z,200);
set(ch,'edgecolor','none');
xlabel('log_{10} (k \tau_1)','FontSize',lab_fontsize);
ylabel('log_{10} (k \tau_2)','FontSize',lab_fontsize);
%zlabel('(\tau_1 \tau_2)^{1/2} <\theta^V (k \tau_1) \theta^V (k \tau_2)> ','FontSize',lab_fontsize);
title(titlename,'FontSize',title_fontsize);
%axis([-2,2.3,-2,2.3])
%axis([2,2.4,2,2.4])
 set(gca,'FontSize',axes_fontsize); hold on;
colorbar;
 set(gca,'FontSize',axes_fontsize); hold on;
%caxis([-0.3 1.0]);
%caxis([-7.3122 47.001]);
%caxis([-0.9095 2.5439]);
%caxis([-0.2788 0.8484]);
%caxis([-0.0229 0.8535]);
%caxis([-4.2977 2.5359]);
caxis

set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperPosition',[ 0 0 6 5]);
outfile = strcat('plots/',filename,'_uetc_',type,'.eps');
print('-depsc2',outfile);



