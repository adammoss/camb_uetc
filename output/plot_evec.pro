pro plot_evec

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
output = 1

filename = ['data/test4']
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

nf = n_elements(filename)

!P.CHARSIZE=1.0
!P.CHARTHICK=3.0
!P.THICK=4
!P.multi=[0,4,4]

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

file = filename[0]+'_uetc_ktau.dat'
nrows = file_lines(file)
ktau = fltarr(1,nrows)
openr,lun,file,/get_lun
readf,lun,ktau
free_lun,lun

if output eq 1 then begin
   set_plot,'ps'
   device,filename='plots/evec.eps'
   device,/color,bits=8,xsize=30,ysize=20	
endif else begin
   set_plot, 'x'
   window,0
endelse

ytitle=textoidl('u')
xtitle=textoidl('k \tau')

file = filename[0]+'_uetc_1_evec.dat'
nrows = file_lines(file)
data1 = fltarr(nrows,nrows)
openr,lun,file,/get_lun
readf,lun,data1
free_lun,lun

file = filename[0]+'_uetc_2_evec.dat'
nrows = file_lines(file)
data2 = fltarr(nrows,nrows)
openr,lun,file,/get_lun
readf,lun,data2
free_lun,lun

file = filename[0]+'_uetc_3_evec.dat'
nrows = file_lines(file)
data3 = fltarr(nrows,nrows)
openr,lun,file,/get_lun
readf,lun,data3
free_lun,lun

;plot,ktau,data1[*,23],xtitle=xtitle,ytitle=ytitle,xthick=6,ythick=6,xstyle=1,ys;tyle=1,/xlog,yrange=[-0.3,0.3],xrange=[1,200]
;oplot,ktau,data2[*,23],color=100
;oplot,ktau,data3[*,23],color=110

for i=0,63 do begin
   plot,ktau,data1[*,i],xtitle=xtitle,ytitle=ytitle,xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=[-0.3,0.3]
   oplot,ktau,data2[*,i],color=100,linestyle=1
   oplot,ktau,data3[*,i],color=110,linestyle=2
endfor

if output eq 1 then begin
   device,/close
endif

end
