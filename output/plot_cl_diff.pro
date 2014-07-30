pro plot_cl_diff

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
output = 1
fid = 'data/default_r3_wmap7_n256_m256_w025'
;filename = ['data/default_r2_n256_m16_w025','data/default_r2_n256_m32_w025','data/default_r2_n256_m64_w025','data/default_r2_n256_m128_w025']
filename = ['data/default_r3_wmap7_n256_m16_w025','data/default_r3_wmap7_n256_m32_w025','data/default_r3_wmap7_n256_m64_w025','data/default_r3_wmap7_n256_m128_w025']
l_range = [2,3000]
yrange = [0.5,1.1]
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

nf = n_elements(filename)

!P.CHARSIZE=1.9
!P.CHARTHICK=3.0
!P.THICK=4
!P.multi=[0,4,3]

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
   device,filename='plots/cls_diff.eps'
   device,/color,bits=8,xsize=30,ysize=20	
  endif else begin
   set_plot, 'x'
   window,0
  endelse

  ;=== Scalars

  ytitle=textoidl('C_l^{TT}/C_l^{TT, fid}')
  xtitle=textoidl('l')	

  file = fid+'_scalCls.dat'
  nrows = file_lines(file)
  data_fid = fltarr(4,nrows)
  openr,lun,file,/get_lun
  readf,lun,data_fid
  free_lun,lun
  data_fid[1,*] = data_fid[1,*]
  data_fid[2,*] = data_fid[2,*]
  data_fid[3,*] = data_fid[3,*]
	
  plot,data_fid[0,*],0*data_fid[1,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=yrange

  for i=0,nf-1 do begin
     file = filename[i]+'_scalCls.dat'
     nrows = file_lines(file)
     data = fltarr(4,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     data[1,*] = data[1,*]
     data[2,*] = data[2,*]
     data[3,*] = data[3,*]
     oplot,data[0,*],data[1,*]/data_fid[1,*],linestyle=i+1,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{EE}/C_l^{EE, fid}')
  xtitle=textoidl('l')		
  plot,data_fid[0,*],data_fid[1,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=yrange

  for i=0,nf-1 do begin
     file = filename[i]+'_scalCls.dat'
     nrows = file_lines(file)
     data = fltarr(4,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     data[1,*] = data[1,*]
     data[2,*] = data[2,*]
     data[3,*] = data[3,*]
     oplot,data[0,*],data[2,*]/data_fid[2,*],linestyle=i+1,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{BB}/C_l^{BB, fid}')
  xtitle=textoidl('l')		
  plot,data_fid[0,*],data_fid[2,*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/nodata,/xlog,yrange=yrange

  ytitle=textoidl('C_l^{TE}/C_l^{TE, fid}')
  xtitle=textoidl('l')		
  plot,data_fid[0,*],data_fid[1,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=yrange

    for i=0,nf-1 do begin
     file = filename[i]+'_scalCls.dat'
     nrows = file_lines(file)
     data = fltarr(4,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     data[1,*] = data[1,*]
     data[2,*] = data[2,*]
     data[3,*] = data[3,*]
     oplot,data[0,*],data[3,*]/data_fid[3,*],linestyle=i+1,color=100+i*10
  endfor

  ;=== Vectors

  file = fid+'_vecCls.dat'
  nrows = file_lines(file)
  data_fid = fltarr(5,nrows)
  openr,lun,file,/get_lun
  readf,lun,data_fid
  free_lun,lun
  data_fid[1,*] = data_fid[1,*]
  data_fid[2,*] = data_fid[2,*]
  data_fid[3,*] = data_fid[3,*]
  data_fid[4,*] = data_fid[4,*]

  ytitle=textoidl('C_l^{TT}/C_l^{TT, fid}')
  xtitle=textoidl('l')		
  plot,data_fid[0,*],data_fid[2,*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=yrange

  for i=0,nf-1 do begin
     file = filename[i]+'_vecCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     data[1,*] = data[1,*]
     data[2,*] = data[2,*]
     data[3,*] = data[3,*]
     data[4,*] = data[4,*]
     oplot,data[0,*],data[1,*]/data_fid[1,*],linestyle=i+1,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{EE}/C_l^{EE, fid}')
  xtitle=textoidl('l')		
  plot,data_fid[0,*],data_fid[2,*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=yrange
 
 for i=0,nf-1 do begin
     file = filename[i]+'_vecCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     data[1,*] = data[1,*]
     data[2,*] = data[2,*]
     data[3,*] = data[3,*]
     data[4,*] = data[4,*]
     oplot,data[0,*],data[2,*]/data_fid[2,*],linestyle=i+1,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{BB}/C_l^{BB, fid}')
  xtitle=textoidl('l')		
  plot,data_fid[0,*],data_fid[1,*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=yrange

 for i=0,nf-1 do begin
     file = filename[i]+'_vecCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     data[1,*] = data[1,*]
     data[2,*] = data[2,*]
     data[3,*] = data[3,*]
     data[4,*] = data[4,*]
     oplot,data[0,*],data[3,*]/data_fid[3,*],linestyle=i+1,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{TE}/C_l^{TE, fid}')
  xtitle=textoidl('l')		
  plot,data_fid[0,*],data_fid[2,*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=yrange
 
 for i=0,nf-1 do begin
     file = filename[i]+'_vecCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     data[1,*] = data[1,*]
     data[2,*] = data[2,*]
     data[3,*] = data[3,*]
     data[4,*] = data[4,*]
     oplot,data[0,*],data[4,*]/data_fid[4,*],linestyle=i+1,color=100+i*10
  endfor

  ;=== Tensors

  file = fid+'_tensCls.dat'
  nrows = file_lines(file)
  data_fid = fltarr(5,nrows)
  openr,lun,file,/get_lun
  readf,lun,data_fid
  free_lun,lun
  data_fid[1,*] = data_fid[1,*]
  data_fid[2,*] = data_fid[2,*]
  data_fid[3,*] = data_fid[3,*]
  data_fid[4,*] = data_fid[4,*]

  ytitle=textoidl('C_l^{TT}/C_l^{TT, fid}')
  xtitle=textoidl('l')		
  plot,data_fid[0,*],data_fid[3,*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=yrange
 
 for i=0,nf-1 do begin
     file = filename[i]+'_tensCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     data[1,*] = data[1,*]
     data[2,*] = data[2,*]
     data[3,*] = data[3,*]
     data[4,*] = data[4,*]
     oplot,data[0,*],data[1,*]/data_fid[1,*],linestyle=i+1,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{EE}/C_l^{EE, fid}')
  xtitle=textoidl('l')		
  plot,data_fid[0,*],data_fid[3,*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=yrange
 
 for i=0,nf-1 do begin
     file = filename[i]+'_tensCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     data[1,*] = data[1,*]
     data[2,*] = data[2,*]
     data[3,*] = data[3,*]
     data[4,*] = data[4,*]
     oplot,data[0,*],data[2,*]/data_fid[2,*],linestyle=i+1,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{BB}/C_l^{BB, fid}')
  xtitle=textoidl('l')		
  plot,data_fid[0,*],data_fid[2,*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=yrange
 
 for i=0,nf-1 do begin
     file = filename[i]+'_tensCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     data[1,*] = data[1,*]
     data[2,*] = data[2,*]
     data[3,*] = data[3,*]
     data[4,*] = data[4,*]
     oplot,data[0,*],data[3,*]/data_fid[3,*],linestyle=i+1,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{TE}/C_l^{TE, fid}')
  xtitle=textoidl('l')		
  plot,data_fid[0,*],data_fid[3,*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=yrange
 
 for i=0,nf-1 do begin
     file = filename[i]+'_tensCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     data[1,*] = data[1,*]
     data[2,*] = data[2,*]
     data[3,*] = data[3,*]
     data[4,*] = data[4,*]
     oplot,data[0,*],data[4,*]/data_fid[4,*],linestyle=i+1,color=100+i*10
  endfor

  if output eq 1 then begin
   device,/close
  endif

  ;exit

end
