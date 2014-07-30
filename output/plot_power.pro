pro plot_power

gmu = 2.0d-7
output = 1
filename = ['data/test']
nf = n_elements(filename)

!P.CHARSIZE=2.0
!P.CHARTHICK=3.0
!P.THICK=4
!P.multi=[0,1,1]

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
   device,filename='plots/power.eps'
   device,/color,bits=8,xsize=20,ysize=18	
  endif else begin
   set_plot, 'x'
   window,0
  endelse

  file = filename[0]+'_uetc_power.dat'
  nrows = file_lines(file)
  data = fltarr(3,nrows)
  openr,lun,file,/get_lun
  readf,lun,data
  free_lun,lun
  data[1,*] = data[1,*]*gmu^2
 
  ytitle=textoidl('Total power [\mu K^2]')
  xtitle=textoidl('Number of eigenmodes')		
  plot,data[0,*],data[1,*]*0.0,xtitle=xtitle,ytitle=ytitle,xrange=[1,100],xthick=6,ythick=6,xstyle=1,ystyle=1,/xlog,yrange=[0.1,30]*10.0^(-10),/ylog

 for i=0,nf-1 do begin
    file = filename[i]+'_uetc_power.dat'
    nrows = file_lines(file)
    data = fltarr(3,nrows)
    openr,lun,file,/get_lun
    readf,lun,data
    free_lun,lun
    data[1,*] = data[1,*]*gmu^2
    data[2,*] = data[2,*]*gmu^2
    data[3,*] = data[3,*]*gmu^2
    oplot,data[0,*],data[1,*],linestyle=i,color=90+i*10
    oplot,data[0,*],data[2,*],linestyle=i,color=90+i*10
    oplot,data[0,*],data[3,*],linestyle=i,color=90+i*10
  endfor

  if output eq 1 then begin
   device,/close
  endif

  ;exit

end
