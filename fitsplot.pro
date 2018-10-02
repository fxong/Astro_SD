;Purpose: quickly plot the derived fitsfile
;Usage: fitsplot, fitsfile, confile=confile, scale=[n1, n2(, n3, n4, n5, n6)], smooth=n, color=n, $
;       position=[n1, n2, n3, n4], /nan, /nofile
;Input:
;  fitsfile: the original fitsfile WITHOUT suffix, in the RGB mode,
;              the input should be [redfits, greenfits, bluefits]
;Input keyword:
;  confile: set a new fitsfile used to plot contours, the default is 
;             to use the original fitsfile to plot contours
;  scale: set the minimun and maximun values of data when plotting, in the RGB mode, 
;           scale should be a six-element array, the default is defined in the procedure
;  smooth: set the smooth degree when plotting, the default is 0
;  color: set the color table when plotting, the default is 0(B-W)
;  position: set the position of figure when plotting, the default is [0.12,0.12,0.9,0.9]
;  nan: set to plot the NAN values in the desired color
;  nofile: set NOT to output the ps file when plotting
;Caution:
;  the default setting is only for quick plotting, the parameters
;    in the flagged area are needed to be set by hand
;  for C180 fitsfile, 'smooth' procedure should be put before 'scale' procedure by hand
;  only the linear scale is available, for other kinds of scale, AICer is better
;History:
;  Completed. 2 Dec 2015.
;  Caution updated and Contour procedure updated. 15 Dec 2015.
;  Add cgplots procedure. 25 Dec 2015.
;  Update the format of colorbar. 11 Jan 2016.
;  Add new keywords, 'position' and 'nofile'. 16 Jan 2016.
;  Scale procedure updated and Bugs fixed. 26 Sep 2016.
;  Set the scales of the X and Y axes to be equal, and remove the NAN procedure when 
;    plotting the single fitsfile under idl85. 29 March 2017.
;Copyright: written by fxong@PMO

pro fitsplot, fitsfile, confile=confile, scale=scale, smooth=smooth, color=color, position=position, nan=nan, nofile=nofile

if n_params() lt 1 then begin
  print,'Syntax - FITSPLOT, fitsfile, [confile= , scale= , smooth= , color= , position= , /nan, /nofile]'
  return
endif
if n_elements(fitsfile) eq 3 then begin
  if ~file_test(fitsfile[0]+'.fits') then begin
    print,'Error: Red file does not exist! Input red file without suffix.'
    return
  endif
  if ~file_test(fitsfile[1]+'.fits') then begin
    print,'Error: Green file does not exist! Input green file without suffix.'
    return
  endif
  if ~file_test(fitsfile[2]+'.fits') then begin
    print,'Error: Blue file does not exist! Input blue file without suffix.'
    return
  endif
endif else begin
  if ~file_test(fitsfile+'.fits') then begin
    print,'Error: Fits file does not exist! Input fits file without suffix.'
    return
  endif
endelse
if keyword_set(confile) then begin
  if ~file_test(confile+'.fits') then begin
    print,'Error: Contour file does not exist! Input contour file without suffix.'
    return
  endif
endif 

if n_elements(fitsfile) eq 3 then begin
  fits_read,fitsfile[0]+'.fits',datred,hdrred
  fits_read,fitsfile[1]+'.fits',datgre,hdrgre
  fits_read,fitsfile[2]+'.fits',datblu,hdrblu
  hdr=hdrred
endif else begin
  fits_read,fitsfile+'.fits',dat,hdr
endelse
if keyword_set(confile) then begin
  fits_read,confile+'.fits',datc,hdrc
endif else begin
  if n_elements(fitsfile) eq 1 then datc=dat
endelse

;axis
xtit=sxpar(hdr,'CTYPE1') & ytit=sxpar(hdr,'CTYPE2')
xmin=sxpar(hdr,'CRVAL1')+(1-sxpar(hdr,'CRPIX1'))*sxpar(hdr,'CDELT1')
xmax=sxpar(hdr,'CRVAL1')+(sxpar(hdr,'NAXIS1')-sxpar(hdr,'CRPIX1'))*sxpar(hdr,'CDELT1')
ymin=sxpar(hdr,'CRVAL2')+(1-sxpar(hdr,'CRPIX2'))*sxpar(hdr,'CDELT2')
ymax=sxpar(hdr,'CRVAL2')+(sxpar(hdr,'NAXIS2')-sxpar(hdr,'CRPIX2'))*sxpar(hdr,'CDELT2')
if keyword_set(confile) then begin
  xminc=sxpar(hdrc,'CRVAL1')+(1-sxpar(hdrc,'CRPIX1'))*sxpar(hdrc,'CDELT1')
  xmaxc=sxpar(hdrc,'CRVAL1')+(sxpar(hdrc,'NAXIS1')-sxpar(hdrc,'CRPIX1'))*sxpar(hdrc,'CDELT1')
  yminc=sxpar(hdrc,'CRVAL2')+(1-sxpar(hdrc,'CRPIX2'))*sxpar(hdrc,'CDELT2')
  ymaxc=sxpar(hdrc,'CRVAL2')+(sxpar(hdrc,'NAXIS2')-sxpar(hdrc,'CRPIX2'))*sxpar(hdrc,'CDELT2')
endif else begin
  if n_elements(fitsfile) eq 1 then xminc=xmin & xmaxc=xmax & yminc=ymin & ymaxc=ymax
endelse

;scale
if keyword_set(scale) then begin
  if n_elements(fitsfile) eq 3 then begin
    mindatr=scale[0] & maxdatr=scale[1]
    datred[where(datred lt mindatr)]=mindatr
    datred[where(datred gt maxdatr)]=maxdatr
    mindatg=scale[2] & maxdatg=scale[3]
    datgre[where(datgre lt mindatg)]=mindatg
    datgre[where(datgre gt maxdatg)]=maxdatg
    mindatb=scale[4] & maxdatb=scale[5]
    datblu[where(datblu lt mindatb)]=mindatb
    datblu[where(datblu gt maxdatb)]=maxdatb
  endif else begin
    mindat=scale[0] & maxdat=scale[1]
    dat[where(dat lt mindat)]=mindat
    dat[where(dat gt maxdat)]=maxdat
  endelse
endif else begin
  if n_elements(fitsfile) eq 3 then begin
    mindatr=ceil(min(datred,/nan))
    maxdatr=floor(max(datred,/nan))
    datred[where(datred lt mindatr)]=mindatr
    datred[where(datred gt maxdatr)]=maxdatr    
    mindatg=ceil(min(datgre,/nan))
    maxdatg=floor(max(datgre,/nan))
    datgre[where(datgre lt mindatg)]=mindatg
    datgre[where(datgre gt maxdatg)]=maxdatg
    mindatb=ceil(min(datblu,/nan))
    maxdatb=floor(max(datblu,/nan))
    datblu[where(datblu lt mindatb)]=mindatb
    datblu[where(datblu gt maxdatb)]=maxdatb
  endif else begin
    mindat=ceil(min(dat,/nan))
    maxdat=floor(max(dat,/nan))
    dat[where(dat lt mindat)]=mindat
    dat[where(dat gt maxdat)]=maxdat
  endelse
endelse

;smooth
if keyword_set(smooth) then begin
  if n_elements(fitsfile) eq 3 then begin
    datred=smooth(datred,smooth,/nan)
    datgre=smooth(datgre,smooth,/nan)
    datblu=smooth(datblu,smooth,/nan)
  endif else begin
    dat=smooth(dat,smooth,/nan)
  endelse
  if keyword_set(confile) then begin
    datc=smooth(datc,smooth,/nan)
  endif else begin
    if n_elements(fitsfile) eq 1 then datc=smooth(datc,smooth,/nan) 
  endelse
endif

;NAN values
if keyword_set(nan) then begin
  if n_elements(fitsfile) eq 3 then begin
    datred=reform(datred,1,sxpar(hdrred,'NAXIS1'),sxpar(hdrred,'NAXIS2'),/overwrite)
    datgre=reform(datgre,1,sxpar(hdrgre,'NAXIS1'),sxpar(hdrgre,'NAXIS2'),/overwrite)
    datblu=reform(datblu,1,sxpar(hdrblu,'NAXIS1'),sxpar(hdrblu,'NAXIS2'),/overwrite)
    datrgb=[datred,datgre,datblu] 
    nanindex=where(finite(datrgb,/nan))
  endif else begin
    ;nanindex=where(finite(dat,/nan)) ;only test on idl85, be cautious when usig other versions of idl!
  endelse
endif else begin
  if n_elements(fitsfile) eq 3 then begin
    datred=reform(datred,1,sxpar(hdrred,'NAXIS1'),sxpar(hdrred,'NAXIS2'),/overwrite)
    datgre=reform(datgre,1,sxpar(hdrgre,'NAXIS1'),sxpar(hdrgre,'NAXIS2'),/overwrite)
    datblu=reform(datblu,1,sxpar(hdrblu,'NAXIS1'),sxpar(hdrblu,'NAXIS2'),/overwrite)
  endif
endelse 

;RGB scale
if n_elements(fitsfile) eq 3 then begin
  imgred=bytscl(datred,min=mindatr,max=maxdatr,/nan)
  imggre=bytscl(datgre,min=mindatg,max=maxdatg,/nan)
  imgblu=bytscl(datblu,min=mindatb,max=maxdatb,/nan)
  image=[imgred,imggre,imgblu]
  if keyword_set(nan) then image[nanindex]=255 ;190, 255
endif else begin
  image=bytscl(dat,min=mindat,max=maxdat,/nan)
  ;if keyword_set(nan) then image[nanindex]=0 ;only test on idl85, be cautious when usig other versions of idl!
endelse

;color
if keyword_set(color) then begin
  loadct,color ;6, 13, 20, 33, 34, 39, 54
endif else begin
  loadct,0
endelse

;position
if keyword_set(position) then begin
  pos=position
endif else begin
  pos=[0.12,0.12,0.9,0.9]
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if ~keyword_set(nofile) then begin
  set_plot,'ps'
  device,filename='fitsplot.eps',/encapsulated,decomposed=1 ;,xsize=21,ysize=29.7
endif

cgplot,[0],position=pos,/isotropic,$
       xrange=[xmin,xmax],$ ;xtickinterval=,xminor=,xtickformat='(f0.1)',xthick=1,$ ;xticklen=,$
       yrange=[ymin,ymax],$ ;ytickinterval=,yminor=,ytickformat='(f0.1)',ythick=1,$ ;yticklen=,$
       xtitle='!6'+string(xtit), ytitle='!6'+string(ytit),charsize=1.2,charthick=1
cgimage,image,$ ;missing_index=0,missing_color='white',$ ;white, grey, black
        xrange=[xmin,xmax],yrange=[ymin,ymax],/overplot,stretch=1 ;1(linear), 4(log)
cgplot,[0],position=pos,/isotropic,$
       xrange=[xmin,xmax],xtickformat='(A1)',$ ;xtickinterval=,xminor=,xthick=1,$ ;xticklen=,$
       yrange=[ymin,ymax],ytickformat='(A1)',$ ;ytickinterval=,yminor=,ythick=1,$ ;yticklen=,$
       charsize=1.2,charthick=1,/noerase,axiscolor='white' ;white, black

;contour
;peak=floor(max(datc,/nan)) ;max(datc,/nan)
;level=[1,2,3] ;[peak*0.1,peak*0.3,peak*0.5,peak*0.7,peak*0.9], [peak*0.2,peak*0.4,peak*0.6,peak*0.8]
;cgcontour,datc,levels=level,c_labels=[0],c_linestyle=[0],c_thick=[1],$
;          c_colors=['white'],/onimage ;position=[]

;colorbar
;cgcolorbar,position=[pos[2],pos[1],pos[2]+0.02,pos[3]],/right,/vertical,format='(i0.0)',$
;           range=[mindat,maxdat],charsize=1.2,charthick=1,title='!6' ;,bottom=min(image[where(image gt 0)]) ;use with very caution, it may change the scale of colorbar!

;marks
;cgplot,[,],[,],psym=0,linestyle=0,thick=1,color='black',/overplot
;cgplot,,,psym=0,linestyle=0,thick=1,color='black',/overplot
;cgplot,,,psym=9,thick=1,symsize=1,symcolor='black',/overplot
;cgtext,,,'!6',charsize=1,charthick=1,color='black' ;,/normal

if ~keyword_set(nofile) then begin
  device,/close_file
  set_plot,'x'
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end

;multiply fitsplot (multiplot.pro)
;Purpose: quickly plot several fitsfiles in one figure
;Copyright: written by fxong@PMO

;@
;@

set_plot,'ps'
device,filename='multiplot.eps',/encapsulated,decomposed=1,xsize=21,ysize=29.7
!P.Multi=[0,1,n]

;,/nofile
;,/nofile

;cgplots,[,],[,],psym=0,linestyle=0,thick=1,color='black',/normal
;cgplots,,,psym=0,linestyle=0,thick=1,color='black',/normal
;cgplots,,,psym=9,thick=1,symsize=1,symcolor='black',/normal

device,/close_file
set_plot,'x'

end
