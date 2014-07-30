pro plot_ttbb

gmu = 2.0d-7
output = 1
fid = 'data/default_n1024_m128_w025'
filename_v = ['data/v02_n1024_m256_w025','data/v04_n1024_m256_w025','data/v06_n1024_m256_w025','data/v08_n1024_m256_w025']
 filename_alp = ['data/alp1_n1024_m256_w025','data/alp15_n1024_m256_w025','data/alp2_n1024_m256_w025','data/alp25_n1024_m256_w025']
 filename_eps = ['data/xi01_n1024_m256_w025','data/xi02_n1024_m256_w025','data/xi03_n1024_m256_w025','data/xi04_n1024_m256_w025']

nf = n_elements(filename_v)
scale_factor_v = fltarr(nf)

!P.CHARSIZE=2.2
!P.CHARTHICK=5.0
!P.THICK=6
!P.multi=[0,3,2]

black = reform([0,0,0],1,3)
tvlct,black,90
red = reform([255,0,0],1,3)
tvlct,red,100
blue = reform([0,0,255],1,3)
tvlct,blue,110
green = reform([0,255,0],1,3)
tvlct,green,120
magenta = reform([255,0,255],1,3)
tvlct,magenta,130
cyan = reform([0,255,255],1,3)
tvlct,cyan,140

 if output eq 1 then begin
   set_plot,'ps'
   device,filename='plots/cls_ttbb.eps'
   device,/color,bits=8,xsize=30,ysize=18	
  endif else begin
   set_plot, 'x'
   window,0
  endelse

  file = fid+'_scalCls.dat'
  nrows = file_lines(file)
  datas = fltarr(4,nrows)
  openr,lun,file,/get_lun
  readf,lun,datas
  free_lun,lun
  datas[1,*] = datas[1,*]*gmu^2
  datas[2,*] = datas[2,*]*gmu^2
  datas[3,*] = datas[3,*]*gmu^2
  file = fid+'_vecCls.dat'
  nrows = file_lines(file)
  datav = fltarr(5,nrows)
  openr,lun,file,/get_lun
  readf,lun,datav
  free_lun,lun
  datav[1,*] = datav[1,*]*gmu^2
  datav[2,*] = datav[2,*]*gmu^2
  datav[3,*] = datav[3,*]*gmu^2
  datav[4,*] = datav[4,*]*gmu^2
  file = fid+'_tensCls.dat'
  nrows = file_lines(file)
  datat = fltarr(5,nrows)
  openr,lun,file,/get_lun
  readf,lun,datat
  free_lun,lun
  datat[1,*] = datat[1,*]*gmu^2
  datat[2,*] = datat[2,*]*gmu^2
  datat[3,*] = datat[3,*]*gmu^2
  datat[4,*] = datat[4,*]*gmu^2
  data_tt = datas[1,*]+datav[1,*]+datat[1,*]

  ytitle=textoidl('C_l^{TT} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  plot,datas[0,*],data_tt[*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=[2,3000],xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=[0.0,550]

 print,'v scaling:'
 for i=0,nf-1 do begin
     file = filename_v[i]+'_scalCls.dat'
     nrows = file_lines(file)
     datas = fltarr(4,nrows)
     openr,lun,file,/get_lun
     readf,lun,datas
     free_lun,lun
     datas[1,*] = datas[1,*]*gmu^2
     datas[2,*] = datas[2,*]*gmu^2
     datas[3,*] = datas[3,*]*gmu^2
     file = filename_v[i]+'_vecCls.dat'
     nrows = file_lines(file)
     datav = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,datav
     free_lun,lun
     datav[1,*] = datav[1,*]*gmu^2
     datav[2,*] = datav[2,*]*gmu^2
     datav[3,*] = datav[3,*]*gmu^2
     datav[4,*] = datav[4,*]*gmu^2
     file = filename_v[i]+'_tensCls.dat'
     nrows = file_lines(file)
     datat = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,datat
     free_lun,lun
     datat[1,*] = datat[1,*]*gmu^2
     datat[2,*] = datat[2,*]*gmu^2
     datat[3,*] = datat[3,*]*gmu^2
     datat[4,*] = datat[4,*]*gmu^2
     scale_factor_v[i] = data_tt[8]/(datas[1,8]+datav[1,8]+datat[1,8])
     oplot,datas[0,*],scale_factor_v[i]*(datas[1,*]+datav[1,*]+datat[1,*]),linestyle=i,color=90+i*10
     print,scale_factor_v[i],sqrt(scale_factor_v[i])
  endfor

 nf = n_elements(filename_alp)
 scale_factor_alp = fltarr(nf)
 
 ytitle=textoidl('C_l^{TT} l(l+1)/(2\pi) [\mu K^2]')
 xtitle=textoidl('l')		
 plot,datas[0,*],data_tt[*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=[2,3000],xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=[0.0,550]

 print,'alpha scaling:'
 for i=0,nf-1 do begin
     file = filename_alp[i]+'_scalCls.dat'
     nrows = file_lines(file)
     datas = fltarr(4,nrows)
     openr,lun,file,/get_lun
     readf,lun,datas
     free_lun,lun
     datas[1,*] = datas[1,*]*gmu^2
     datas[2,*] = datas[2,*]*gmu^2
     datas[3,*] = datas[3,*]*gmu^2
     file = filename_alp[i]+'_vecCls.dat'
     nrows = file_lines(file)
     datav = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,datav
     free_lun,lun
     datav[1,*] = datav[1,*]*gmu^2
     datav[2,*] = datav[2,*]*gmu^2
     datav[3,*] = datav[3,*]*gmu^2
     datav[4,*] = datav[4,*]*gmu^2
     file = filename_alp[i]+'_tensCls.dat'
     nrows = file_lines(file)
     datat = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,datat
     free_lun,lun
     datat[1,*] = datat[1,*]*gmu^2
     datat[2,*] = datat[2,*]*gmu^2
     datat[3,*] = datat[3,*]*gmu^2
     datat[4,*] = datat[4,*]*gmu^2
     scale_factor_alp[i] = data_tt[8]/(datas[1,8]+datav[1,8]+datat[1,8])
     oplot,datas[0,*],scale_factor_alp[i]*(datas[1,*]+datav[1,*]+datat[1,*]),linestyle=i,color=90+i*10
     print,scale_factor_alp[i],sqrt(scale_factor_alp[i])
  endfor

 nf = n_elements(filename_eps)
 scale_factor_eps = fltarr(nf)

 ytitle=textoidl('C_l^{TT} l(l+1)/(2\pi) [\mu K^2]')
 xtitle=textoidl('l')		
 plot,datas[0,*],data_tt[*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=[2,3000],xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=[0.0,550]

 print,'epsilon scaling:'
 for i=0,nf-1 do begin
     file = filename_eps[i]+'_scalCls.dat'
     nrows = file_lines(file)
     datas = fltarr(4,nrows)
     openr,lun,file,/get_lun
     readf,lun,datas
     free_lun,lun
     datas[1,*] = datas[1,*]*gmu^2
     datas[2,*] = datas[2,*]*gmu^2
     datas[3,*] = datas[3,*]*gmu^2
     file = filename_eps[i]+'_vecCls.dat'
     nrows = file_lines(file)
     datav = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,datav
     free_lun,lun
     datav[1,*] = datav[1,*]*gmu^2
     datav[2,*] = datav[2,*]*gmu^2
     datav[3,*] = datav[3,*]*gmu^2
     datav[4,*] = datav[4,*]*gmu^2
     file = filename_eps[i]+'_tensCls.dat'
     nrows = file_lines(file)
     datat = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,datat
     free_lun,lun
     datat[1,*] = datat[1,*]*gmu^2
     datat[2,*] = datat[2,*]*gmu^2
     datat[3,*] = datat[3,*]*gmu^2
     datat[4,*] = datat[4,*]*gmu^2
     scale_factor_eps[i] = data_tt[8]/(datas[1,8]+datav[1,8]+datat[1,8])
     oplot,datas[0,*],scale_factor_eps[i]*(datas[1,*]+datav[1,*]+datat[1,*]),linestyle=i,color=90+i*10
     print,scale_factor_eps[i],sqrt(scale_factor_eps[i])
  endfor
 
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  nf = n_elements(filename_v)

  ytitle=textoidl('C_l^{BB} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  plot,datas[0,*],data_tt[*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=[2,3000],xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=[0.0001,1]

 for i=0,nf-1 do begin
     file = filename_v[i]+'_vecCls.dat'
     nrows = file_lines(file)
     datav = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,datav
     free_lun,lun
     datav[1,*] = datav[1,*]*gmu^2
     datav[2,*] = datav[2,*]*gmu^2
     datav[3,*] = datav[3,*]*gmu^2
     datav[4,*] = datav[4,*]*gmu^2
     file = filename_v[i]+'_tensCls.dat'
     nrows = file_lines(file)
     datat = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,datat
     free_lun,lun
     datat[1,*] = datat[1,*]*gmu^2
     datat[2,*] = datat[2,*]*gmu^2
     datat[3,*] = datat[3,*]*gmu^2
     datat[4,*] = datat[4,*]*gmu^2
     oplot,datav[0,*],scale_factor_v[i]*(datav[3,*]+datat[3,*]),linestyle=i,color=90+i*10
  endfor

 nf = n_elements(filename_alp)

 ytitle=textoidl('C_l^{BB} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  plot,datas[0,*],data_tt[*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=[2,3000],xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=[0.0001,1]

 for i=0,nf-1 do begin
    file = filename_alp[i]+'_vecCls.dat'
     nrows = file_lines(file)
     datav = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,datav
     free_lun,lun
     datav[1,*] = datav[1,*]*gmu^2
     datav[2,*] = datav[2,*]*gmu^2
     datav[3,*] = datav[3,*]*gmu^2
     datav[4,*] = datav[4,*]*gmu^2
     file = filename_alp[i]+'_tensCls.dat'
     nrows = file_lines(file)
     datat = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,datat
     free_lun,lun
     datat[1,*] = datat[1,*]*gmu^2
     datat[2,*] = datat[2,*]*gmu^2
     datat[3,*] = datat[3,*]*gmu^2
     datat[4,*] = datat[4,*]*gmu^2
     oplot,datav[0,*],scale_factor_alp[i]*(datav[3,*]+datat[3,*]),linestyle=i,color=90+i*10
  endfor

 nf = n_elements(filename_eps)

 ytitle=textoidl('C_l^{BB} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  plot,datas[0,*],data_tt[*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=[2,3000],xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=[0.0001,1]

 for i=0,nf-1 do begin
     file = filename_eps[i]+'_vecCls.dat'
     nrows = file_lines(file)
     datav = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,datav
     free_lun,lun
     datav[1,*] = datav[1,*]*gmu^2
     datav[2,*] = datav[2,*]*gmu^2
     datav[3,*] = datav[3,*]*gmu^2
     datav[4,*] = datav[4,*]*gmu^2
     file = filename_eps[i]+'_tensCls.dat'
     nrows = file_lines(file)
     datat = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,datat
     free_lun,lun
     datat[1,*] = datat[1,*]*gmu^2
     datat[2,*] = datat[2,*]*gmu^2
     datat[3,*] = datat[3,*]*gmu^2
     datat[4,*] = datat[4,*]*gmu^2
     oplot,datav[0,*],scale_factor_eps[i]*(datav[3,*]+datat[3,*]),linestyle=i,color=90+i*10
  endfor

  if output eq 1 then begin
   device,/close
  endif

  ;exit

end
