;Purpose: clip the fits into subfits of the desired area, 
;           also provide pixel values of the desired coordinates
;Usage: fitsclip, fitsfile, x=[x1, x2], y=[y1, y2], z=[z1, z2], /onlypix
;Input:
;  fitsfile: the original fitsfile WITHOUT suffix
;Input keyword:
;  x: set the range [x1, x2] to clip the datacube along X axis
;  y: set the range [y1, y2] to clip the datacube along Y axis
;  z: set the range [z1, z2] to clip the datacube along Z axis
;  onlypix: set to only show the pixel values of the desired coordinates, 
;             without clipping the fits
;Caution:
;  the unit of x, y and z should refer to the head of fitsfile
;History:
;  Completed. 22 Nov 2015.
;Copyright: written by fxong@PMO

function coord2pix, hdr, value, dim, ds9pix=ds9pix
value=float(value)
dim=strcompress(string(dim),/remove_all)
pix=(value-sxpar(hdr,'CRVAL'+dim))/sxpar(hdr,'CDELT'+dim)+sxpar(hdr,'CRPIX'+dim)
if keyword_set(ds9pix) then return,pix else return,pix-1
end

pro fitsclip, fitsfile, x=x, y=y, z=z, onlypix=onlypix

if n_params() lt 1 then begin
  print,'Syntax - FITSCLIP, fitsfile, [x= , y= , z= , /onlypix]'
  return
endif
if ~file_test(fitsfile+'.fits') then begin
  print,'Error: Fits file does not exist! Input fits file without suffix.'
  return
endif

fits_read,fitsfile+'.fits',dat,hdr

if keyword_set(x) then begin
  x1=round(coord2pix(hdr,x[0],1))
  x2=round(coord2pix(hdr,x[1],1))
  num1=x2-x1+1
  sxaddpar,hdr,'CRVAL1',x[0]
  sxaddpar,hdr,'CRPIX1',1
endif else begin
  x1=0
  num1=sxpar(hdr,'NAXIS1')
endelse

if keyword_set(y) then begin
  y1=round(coord2pix(hdr,y[0],2))
  y2=round(coord2pix(hdr,y[1],2))
  num2=y2-y1+1
  sxaddpar,hdr,'CRVAL2',y[0]
  sxaddpar,hdr,'CRPIX2',1
endif else begin
  y1=0
  num2=sxpar(hdr,'NAXIS2')
endelse

if sxpar(hdr,'NAXIS') eq 3 then begin
  if keyword_set(z) then begin
    z1=round(coord2pix(hdr,z[0],3))
    z2=round(coord2pix(hdr,z[1],3))
    num3=z2-z1+1
    sxaddpar,hdr,'CRVAL3',z[0]
    sxaddpar,hdr,'CRPIX3',1
  endif else begin
    z1=0
    num3=sxpar(hdr,'NAXIS3')
  endelse
endif

sxdelpar,hdr,['NAXIS','CTYPE','CRVAL','CDELT','CRPIX','CROTA']+'4' ;MWISP fits
sxdelpar,hdr,['ORIGIN','DATE','RMS','HISTORY']                     ;MWISP fits
sxaddhist,'Clip file: '+fitsfile+'.fits',hdr,/comment

if keyword_set(onlypix) then begin

  if keyword_set(x) then begin
    print,'The pixel value of',x[0],'  on X axis is',x1,'  (in ds9 is',x1-1,')'
    print,'The pixel value of',x[1],'  on X axis is',x2,'  (in ds9 is',x2-1,')'
  endif
  if keyword_set(y) then begin
    print,'The pixel value of',y[0],'  on Y axis is',y1,'  (in ds9 is',y1-1,')'
    print,'The pixel value of',y[1],'  on Y axis is',y2,'  (in ds9 is',y2-1,')'
  endif
  if keyword_set(z) then begin
    print,'The pixel value of',z[0],'  on Z axis is',z1,'  (in ds9 is',z1-1,')'
    print,'The pixel value of',z[1],'  on Z axis is',z2,'  (in ds9 is',z2-1,')'
  endif

endif else begin

  if sxpar(hdr,'NAXIS') eq 3 then begin
    sldat=dblarr(num1,num2,num3)
    for i=0,num1-1 do begin
      for j=0,num2-1 do begin
        for k=0,num3-1 do begin
          sldat[i,j,k]=dat[i+x1,j+y1,k+z1]
        endfor
      endfor
    endfor
  endif else begin
    sldat=dblarr(num1,num2)
    for i=0,num1-1 do begin
      for j=0,num2-1 do begin
        sldat[i,j]=dat[i+x1,j+y1]
      endfor
    endfor
  endelse
  fits_write,fitsfile+'_clip.fits',sldat,hdr

  if keyword_set(x) then begin
    print,'The pixel value of',x[0],'  on X axis is',x1,'  (in ds9 is',x1-1,')'
    print,'The pixel value of',x[1],'  on X axis is',x2,'  (in ds9 is',x2-1,')'
  endif
  if keyword_set(y) then begin
    print,'The pixel value of',y[0],'  on Y axis is',y1,'  (in ds9 is',y1-1,')'
    print,'The pixel value of',y[1],'  on Y axis is',y2,'  (in ds9 is',y2-1,')'
  endif
  if keyword_set(z) then begin
    print,'The pixel value of',z[0],'  on Z axis is',z1,'  (in ds9 is',z1-1,')'
    print,'The pixel value of',z[1],'  on Z axis is',z2,'  (in ds9 is',z2-1,')'
  endif

endelse

end
