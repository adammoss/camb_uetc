pro plot_vt

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
gmu = 2.0d-7
output = 1

;filename = ['data/xi02_n1024_m256_w025','data/xi04_n1024_m256_w025']
;v1 = 0.65
;eps1 = 0.2
;alp1 = 1.9
;v2 = 0.65
;eps2 = 0.4
;alp2 = 1.9

;filename = ['data/v04_n1024_m256_w025','data/v08_n1024_m256_w025']
;v1 = 0.4
;eps1 = 0.13
;alp1 = 1.9
;v2 = 0.8
;eps2 = 0.13
;alp2 = 1.9

filename = ['data/alp1_n1024_m256_w025','data/alp2_n1024_m256_w025']
v1 = 0.65
eps1 = 0.13
alp1 = 1.0
v2 = 0.65
eps2 = 0.13
alp2 = 2.0

scale1 = (1+v1^2*(alp1^2-2.0)+v1^4*(1.0-alp1^2+alp1^4))/(eps1*alp1^2*(1-v1^2))
scale2 = (1+v2^2*(alp2^2-2.0)+v2^4*(1.0-alp2^2+alp2^4))/(eps2*alp2^2*(1-v2^2))

normalization = [gmu,gmu]

l_range = [2,3000]
tt_s_range = [0,300]
ee_s_range = [0.0001,10]
te_s_range = [-3,0.5]
tt_v_range = [0.0,250]
ee_v_range = [0.00005,0.1]
bb_v_range = [0.0001,1]
te_v_range = [0.0,0.5]
tt_t_range = [0.0,10]
ee_t_range = [0.00001,0.01]
bb_t_range = [0.00001,0.01]
te_t_range = [-0.02,0.015]
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!P.CHARSIZE=1.9
!P.CHARTHICK=3.0
!P.THICK=4
!P.multi=[0,3,2]

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
yellow = reform([255,255,0],1,3)
tvlct,yellow,150

 if output eq 1 then begin
   set_plot,'ps'
   device,filename='plots/cls_vt.eps'
   device,/color,bits=8,xsize=30,ysize=20	
  endif else begin
   set_plot, 'x'
   window,0
  endelse

  ;=== Vectors

  file = filename[0]+'_vecCls.dat'
  nrows = file_lines(file)
  data1 = fltarr(5,nrows)
  openr,lun,file,/get_lun
  readf,lun,data1
  free_lun,lun
  norm = normalization[0]
  data1[1,*] = data1[1,*]*norm^2
  data1[2,*] = data1[2,*]*norm^2
  data1[3,*] = data1[3,*]*norm^2
  data1[4,*] = data1[4,*]*norm^2
  file = filename[1]+'_vecCls.dat'
  nrows = file_lines(file)
  data2 = fltarr(5,nrows)
  openr,lun,file,/get_lun
  readf,lun,data2
  free_lun,lun
  norm = normalization[0]
  data2[1,*] = data2[1,*]*norm^2
  data2[2,*] = data2[2,*]*norm^2
  data2[3,*] = data2[3,*]*norm^2
  data2[4,*] = data2[4,*]*norm^2
  

  ytitle=textoidl('C_l^{TT} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  plot,data1[0,*],data1[1,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=tt_v_range
  oplot,data2[0,*],data2[1,*]*scale1/scale2,color=100,linestyle=1
 
  ytitle=textoidl('C_l^{EE} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  plot,data1[0,*],data1[2,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=ee_v_range
  oplot,data2[0,*],data2[2,*]*scale1/scale2,color=100,linestyle=1

  ytitle=textoidl('C_l^{BB} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  plot,data1[0,*],data1[3,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=bb_v_range
  oplot,data2[0,*],data2[3,*]*scale1/scale2,color=100,linestyle=1


  file = filename[0]+'_tensCls.dat'
  nrows = file_lines(file)
  data1 = fltarr(5,nrows)
  openr,lun,file,/get_lun
  readf,lun,data1
  free_lun,lun
  norm = normalization[0]
  data1[1,*] = data1[1,*]*norm^2
  data1[2,*] = data1[2,*]*norm^2
  data1[3,*] = data1[3,*]*norm^2
  data1[4,*] = data1[4,*]*norm^2
  file = filename[1]+'_tensCls.dat'
  nrows = file_lines(file)
  data2 = fltarr(5,nrows)
  openr,lun,file,/get_lun
  readf,lun,data2
  free_lun,lun
  norm = normalization[0]
  data2[1,*] = data2[1,*]*norm^2
  data2[2,*] = data2[2,*]*norm^2
  data2[3,*] = data2[3,*]*norm^2
  data2[4,*] = data2[4,*]*norm^2
  

  ytitle=textoidl('C_l^{TT} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  plot,data1[0,*],data1[1,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=tt_t_range
  oplot,data2[0,*],data2[1,*]*scale1/scale2,color=100,linestyle=1
 
  ytitle=textoidl('C_l^{EE} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  plot,data1[0,*],data1[2,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=ee_t_range
  oplot,data2[0,*],data2[2,*]*scale1/scale2,color=100,linestyle=1

  ytitle=textoidl('C_l^{BB} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  plot,data1[0,*],data1[3,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=bb_t_range
  oplot,data2[0,*],data2[3,*]*scale1/scale2,color=100,linestyle=1


  if output eq 1 then begin
   device,/close
  endif

  ;exit

end
