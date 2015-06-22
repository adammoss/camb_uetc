pro plot_cl

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
gmu = 2.0d-7
output = 1
plot_cmbact = 0

;cmbact = ['fid_wmap7_500_L099_n500_1','fid_wmap7_500_L099_n500_2','fid_wmap7_500_L099_n500_3','fid_wmap7_500_L099_n500_4']

cmbact = ['test_L099']
;cmbact = []

filename = ['data/test','data/test2']
;filename = ['data/test13','data/test14']
;filename = ['data/v02_r5_wmap7_n1024_m256_w025']

;normalization = [gmu,gmu,gmu,gmu,gmu]
normalization = [gmu,gmu]

l_range = [2,3000]
tt_s_range = [0,100]
ee_s_range = [0.0001,10]
te_s_range = [-3,0.5]
tt_v_range = [0.0,100]
ee_v_range = [0.00005,0.1]
bb_v_range = [0.0001,1]
te_v_range = [0.0,0.5]
tt_t_range = [0.0,6]
ee_t_range = [0.00001,0.02]
bb_t_range = [0.00001,0.01]
te_t_range = [-0.02,0.015]

cmbact_fact = 1.0

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

nf = n_elements(filename)
ncmbact = n_elements(cmbact)
nnorm = n_elements(normalization)
if nf ne nnorm then begin
   print,'Gmu normalization option needs to be same size as input spectrum'
   stop
endif

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
yellow = reform([255,255,0],1,3)
tvlct,yellow,150

file = 'cmbact_data/cl_tt_'+cmbact[0]+'.d'
nrows = file_lines(file)
data = fltarr(4,nrows)
openr,lun,file,/get_lun
readf,lun,data
free_lun,lun

nlmax = n_elements(data[1,*])
data_tt = fltarr(4,nlmax)
data_ee = fltarr(4,nlmax)
data_bb = fltarr(3,nlmax)
data_te = fltarr(4,nlmax)

for i=0,ncmbact-1 do begin
   file = 'cmbact_data/cl_tt_'+cmbact[i]+'.d'
   nrows = file_lines(file)
   data = fltarr(4,nrows)
   openr,lun,file,/get_lun
   readf,lun,data
   free_lun,lun

   data_tt[0,*] += data[0,*]
   data_tt[1,*] += data[1,*]*(gmu/1.1d-6)^2*cmbact_fact
   data_tt[2,*] += data[2,*]*(gmu/1.1d-6)^2*cmbact_fact
   data_tt[3,*] += data[3,*]*(gmu/1.1d-6)^2*cmbact_fact

   file = 'cmbact_data/cl_ee_'+cmbact[i]+'.d'
   nrows = file_lines(file)
   data = fltarr(4,nrows)
   openr,lun,file,/get_lun
   readf,lun,data
   free_lun,lun

   data_ee[0,*] += data[0,*]
   data_ee[1,*] += data[1,*]*(gmu/1.1d-6)^2*cmbact_fact
   data_ee[2,*] += data[2,*]*(gmu/1.1d-6)^2*cmbact_fact
   data_ee[3,*] += data[3,*]*(gmu/1.1d-6)^2*cmbact_fact

   file = 'cmbact_data/cl_bb_'+cmbact[i]+'.d'
   nrows = file_lines(file)
   data = fltarr(3,nrows)
   openr,lun,file,/get_lun
   readf,lun,data
   free_lun,lun
   
   data_bb[0,*] += data[0,*]
   data_bb[1,*] += data[1,*]*(gmu/1.1d-6)^2*cmbact_fact
   data_bb[2,*] += data[2,*]*(gmu/1.1d-6)^2*cmbact_fact
   
   file = 'cmbact_data/cl_te_'+cmbact[i]+'.d'
   nrows = file_lines(file)
   data = fltarr(4,nrows)
   openr,lun,file,/get_lun
   readf,lun,data
   free_lun,lun

   data_te[0,*] += data[0,*]
   data_te[1,*] += data[1,*]*(gmu/1.1d-6)^2*cmbact_fact
   data_te[2,*] += -data[2,*]*(gmu/1.1d-6)^2*cmbact_fact
   data_te[3,*] += data[3,*]*(gmu/1.1d-6)^2*cmbact_fact

end

data_tt = data_tt/ncmbact
data_ee = data_ee/ncmbact
data_bb = data_bb/ncmbact
data_te = data_te/ncmbact

print,'Correcting for minus sign in CMBACT vector TE spectra'

 if output eq 1 then begin
   set_plot,'ps'
   device,filename='plots/cls.eps'
   device,/color,bits=8,xsize=30,ysize=20	
  endif else begin
   set_plot, 'x'
   window,0
  endelse

  ;=== Scalars

  ytitle=textoidl('C_l^{TT} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  if plot_cmbact eq 1 then begin
     plot,data_tt[0,*],data_tt[1,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=tt_s_range
  endif else begin
     plot,data_tt[0,*],data_tt[1,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=tt_s_range,/nodata
  endelse

  for i=0,nf-1 do begin
     file = filename[i]+'_scalCls.dat'
     nrows = file_lines(file)
     data = fltarr(4,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     norm = normalization[i]
     data[1,*] = data[1,*]*norm^2
     data[2,*] = data[2,*]*norm^2
     data[3,*] = data[3,*]*norm^2
     oplot,data[0,*],data[1,*],linestyle=0,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{EE} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')	
  if plot_cmbact eq 1 then begin
     plot,data_ee[0,*],data_ee[1,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=ee_s_range
  endif else begin
     plot,data_ee[0,*],data_ee[1,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=ee_s_range,/nodata
  endelse

  for i=0,nf-1 do begin
     file = filename[i]+'_scalCls.dat'
     nrows = file_lines(file)
     data = fltarr(4,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     norm = normalization[i]
     data[1,*] = data[1,*]*norm^2
     data[2,*] = data[2,*]*norm^2
     data[3,*] = data[3,*]*norm^2
     oplot,data[0,*],data[2,*],linestyle=0,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{BB} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  plot,data_bb[0,*],data_bb[2,*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/nodata,/xlog

  ytitle=textoidl('C_l^{TE} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')	
  if plot_cmbact eq 1 then begin
     plot,data_te[0,*],data_te[1,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=te_s_range
  endif else begin
     plot,data_te[0,*],data_te[1,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=te_s_range,/nodata
  endelse

    for i=0,nf-1 do begin
     file = filename[i]+'_scalCls.dat'
     nrows = file_lines(file)
     data = fltarr(4,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     norm = normalization[i]
     data[1,*] = data[1,*]*norm^2
     data[2,*] = data[2,*]*norm^2
     data[3,*] = data[3,*]*norm^2
     oplot,data[0,*],data[3,*],linestyle=0,color=100+i*10
  endfor

  ;=== Vectors

  ytitle=textoidl('C_l^{TT} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')	
  if plot_cmbact eq 1 then begin
     plot,data_tt[0,*],data_tt[2,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=tt_v_range
  endif else begin
     plot,data_tt[0,*],data_tt[2,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=tt_v_range,/nodata
  endelse

  for i=0,nf-1 do begin
     file = filename[i]+'_vecCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     norm = normalization[i]
     data[1,*] = data[1,*]*norm^2
     data[2,*] = data[2,*]*norm^2
     data[3,*] = data[3,*]*norm^2
     data[4,*] = data[4,*]*norm^2
     oplot,data[0,*],data[1,*],linestyle=0,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{EE} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')	
  if plot_cmbact eq 1 then begin
     plot,data_ee[0,*],data_ee[2,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=ee_v_range
  endif else begin
     plot,data_ee[0,*],data_ee[2,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=ee_v_range,/nodata
  endelse

 for i=0,nf-1 do begin
     file = filename[i]+'_vecCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     norm = normalization[i]
     data[1,*] = data[1,*]*norm^2
     data[2,*] = data[2,*]*norm^2
     data[3,*] = data[3,*]*norm^2
     data[4,*] = data[4,*]*norm^2
     oplot,data[0,*],data[2,*],linestyle=0,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{BB} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')	
  if plot_cmbact eq 1 then begin
     plot,data_bb[0,*],data_bb[1,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=bb_v_range
  endif else begin
     plot,data_bb[0,*],data_bb[1,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=bb_v_range,/nodata
  endelse

 for i=0,nf-1 do begin
     file = filename[i]+'_vecCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     norm = normalization[i]
     data[1,*] = data[1,*]*norm^2
     data[2,*] = data[2,*]*norm^2
     data[3,*] = data[3,*]*norm^2
     data[4,*] = data[4,*]*norm^2
     oplot,data[0,*],data[3,*],linestyle=0,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{TE} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  if plot_cmbact eq 1 then begin
     plot,data_te[0,*],data_te[2,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=te_v_range
  endif else begin
     plot,data_te[0,*],data_te[2,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=te_v_range,/nodata
  endelse

 for i=0,nf-1 do begin
     file = filename[i]+'_vecCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     norm = normalization[i]
     data[1,*] = data[1,*]*norm^2
     data[2,*] = data[2,*]*norm^2
     data[3,*] = data[3,*]*norm^2
     data[4,*] = data[4,*]*norm^2
     oplot,data[0,*],data[4,*],linestyle=0,color=100+i*10
  endfor

  ;=== Tensors

  ytitle=textoidl('C_l^{TT} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  if plot_cmbact eq 1 then begin
     plot,data_tt[0,*],data_tt[3,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=tt_t_range
  endif else begin
     plot,data_tt[0,*],data_tt[3,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=tt_t_range,/nodata
  endelse

 for i=0,nf-1 do begin
     file = filename[i]+'_tensCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     norm = normalization[i]
     data[1,*] = data[1,*]*norm^2
     data[2,*] = data[2,*]*norm^2
     data[3,*] = data[3,*]*norm^2
     data[4,*] = data[4,*]*norm^2
     oplot,data[0,*],data[1,*],linestyle=0,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{EE} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  if plot_cmbact eq 1 then begin
     plot,data_ee[0,*],data_ee[3,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=ee_t_range
  endif else begin
     plot,data_ee[0,*],data_ee[3,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=ee_t_range,/nodata
  endelse

 for i=0,nf-1 do begin
     file = filename[i]+'_tensCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     norm = normalization[i]
     data[1,*] = data[1,*]*norm^2
     data[2,*] = data[2,*]*norm^2
     data[3,*] = data[3,*]*norm^2
     data[4,*] = data[4,*]*norm^2
     oplot,data[0,*],data[2,*],linestyle=0,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{BB} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  if plot_cmbact eq 1 then begin
     plot,data_bb[0,*],data_bb[2,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=bb_t_range
  endif else begin
     plot,data_bb[0,*],data_bb[2,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,/ylog,yrange=bb_t_range,/nodata
  endelse

 for i=0,nf-1 do begin
     file = filename[i]+'_tensCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     norm = normalization[i]
     data[1,*] = data[1,*]*norm^2
     data[2,*] = data[2,*]*norm^2
     data[3,*] = data[3,*]*norm^2
     data[4,*] = data[4,*]*norm^2
     oplot,data[0,*],data[3,*],linestyle=0,color=100+i*10
  endfor

  ytitle=textoidl('C_l^{TE} l(l+1)/(2\pi) [\mu K^2]')
  xtitle=textoidl('l')		
  if plot_cmbact eq 1 then begin
     plot,data_te[0,*],data_te[3,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=te_t_range
  endif else begin
     plot,data_te[0,*],data_te[3,*],xtitle=xtitle,ytitle=ytitle,xrange=l_range,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=te_t_range,/nodata
  endelse

 for i=0,nf-1 do begin
     file = filename[i]+'_tensCls.dat'
     nrows = file_lines(file)
     data = fltarr(5,nrows)
     openr,lun,file,/get_lun
     readf,lun,data
     free_lun,lun
     norm = normalization[i]
     data[1,*] = data[1,*]*norm^2
     data[2,*] = data[2,*]*norm^2
     data[3,*] = data[3,*]*norm^2
     data[4,*] = data[4,*]*norm^2
     oplot,data[0,*],data[4,*],linestyle=0,color=100+i*10
  endfor

  if output eq 1 then begin
   device,/close
  endif

  ;exit

end
