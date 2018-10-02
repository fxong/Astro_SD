;Purpose: Extract position-velocity or position-position map from a datacube
;Usage: ppvslice, fitsfile, [x1, x2, x3,...], [y1, y2, y3,...], v=[v1, v2], $
;          step=n, width=n, /gal, /spline, /track, /pps
;       ppvslice, fitsfile, 'path.cat', v=[v1, v2], step=n, width=n, $
;          /gal, /spline, /track, /pps
;Input:
;  fitsfile: a string indicates the file name of datacube without suffix
;  a:	if 'd' is not present, 'a' is a string indicates the file name of the path, in "x y" format
;  a,d:	coordinate of the path, a=[x1,x2,x3,...], d=[y1,y2,x3,...]
;Input keyword:
;  v:      set the velocity range of p-v map, the default is the original velocity range
;  step:	resample the input slice path with interval of step (in pixel unit), default is 0.5
;  width:	the width (in pixel unit) of the belt to be averaged, default is 1
;  gal:	set to 1 if the given a, d are in galactic system
;  spline:	set to 1 to do spline interpolation to the input path, the path will become smooth spline
;  track:  set to 1 to derive the track info of the path (and the outline of the belt)
;  pps:    set to 1 to derive the position-position map of a datacube
;Caution:
;  only accept fits in celestial and Galactic system with CDELT1=CDELT2
;  the default is to derive the position-velocity map of a 3D datacube, $
;    set 'pps' to 1 to derive the position-position map
;  keyword 'pps' only works for the 2D datacube, use with very caution!!!
;  the 'CDELT' of output fitsfile is half of the original 'CDELT' along the width axis when keyword 'pps' is set
;  the unit of v1 and v2 is km/s
;History:
;  Completed. 02 Dec 2016.
;Copyright: modified by fxong@PMO

function coord2pix, hdr, value, dim, ds9pix=ds9pix
value=float(value)
dim=strcompress(string(dim),/remove_all)
pix=(value-sxpar(hdr,'CRVAL'+dim))/sxpar(hdr,'CDELT'+dim)+sxpar(hdr,'CRPIX'+dim)
if keyword_set(ds9pix) then return,pix else return,pix-1
end

pro clipvel, olddat, oldhdr, newdat, newhdr, vel1, vel2
;clip the datacube with the given velocity range
if sxpar(oldhdr,'NAXIS') eq 3 then begin
  x1=0
  num1=sxpar(oldhdr,'NAXIS1')
  y1=0
  num2=sxpar(oldhdr,'NAXIS2')
  z1=round(coord2pix(oldhdr,vel1,3))
  z2=round(coord2pix(oldhdr,vel2,3))
  num3=z2-z1+1
  newdat=dblarr(num1,num2,num3)
  for i=0,num1-1 do begin
    for j=0,num2-1 do begin
      for k=0,num3-1 do begin
        newdat[i,j,k]=olddat[i+x1,j+y1,k+z1]
      endfor
    endfor
  endfor

  newhdr=oldhdr
  sxaddpar,newhdr,'NAXIS3',num3
  sxaddpar,newhdr,'CRVAL3',vel1
  sxaddpar,newhdr,'CRPIX3',1
endif else begin
  print,'Error: velocity clip only works for 3D datacube! Check the input datacube.'
  newdat=olddat
  newhdr=oldhdr
endelse

end

;pvslice, pvshow
;by ShaoboZhang
;History:
;Apr,26,2012,v1.0
;May,16,2012,v1.1
;  PVSHOW: delete vrange keyword, add keyword: line,transpose,arcmin,log_scale.
;          now can plot l-v or b-v map created by cubemoment.
;  PVSLICE: move the reference pixel of position to the center of the axis
;Jan,04,2016,v2.0
;  PVSLICE: rewrite the procedure, add WIDTH to the slice, correct the path resampling error

pro ppvslice, fitsfile, a, d, v=v, width=width, step=step, gal=gal, spline=spline, track=track, pps=pps;, kernel=kernel
;Extract position-velocity map from a data cube
;Only accept fits in celestial and Galactic system with CDELT1=CDELT2
;Input:
;	fitsfile: a string indicates the file name of datacube
;	a:	if 'd' is not present, 'a' is a string indicates the file name of the path, in "x y" format
;	a,d:	coordinate of the path, a=[x1,x2,x3,...], d=[y1,y2,x3,...]
;Input keyword:
;	width:	the width (in pixel unit) of the belt to be averaged, default is 1
;	step:	resample the input slice path with interval of step (in pixel unit), default is 0.5
;	gal:	set to 1 if the given a, d are in galactic system
;	spline:	set to 1 to do spline interpolation to the input path, the path will become smooth spline
;	kernel: not work for now
;Output:
;	"pvslice.fits": pv map
;	"pvslice.path": the resampled path
;	"pvslice.belt": the outline of the belt
;Usage: pvslice, 'XXX.fits','path.cat',/gal
;	pvslice, 'XXX.fits',[x1,x2,x3],[y1,y2,y3],width=5,/spline

;regular input check
if n_params() le 1 then begin
	print, 'Syntax - PPVSLICE, fitsfile, catalog, [v= , step= , width= , /gal, /spline, /track, /pps]'
	print, 'Syntax - PPVSLICE, fitsfile, a, d, [v= , step= , width= , /gal, /spline, /track, /pps]'
	return
endif
if ~file_test(fitsfile+'.fits') then begin
	print,'Error: fits file does not exist!  Input fits file without suffix.'
	return
endif

;read fits (and clip the velocity range)
if keyword_set(v) then begin
  print,'Run the clipvel'
  fits_read,fitsfile+'.fits',olddat,oldhdr
  factor=abs(sxpar(oldhdr,'CDELT3') lt 100)?(1d):(1000d)
  clipvel, olddat, oldhdr, dat, hdr, v[0]*factor, v[1]*factor  ;clip the velocity range
endif else begin
  fits_read,fitsfile+'.fits',dat,hdr
endelse

;read catalog
if n_params() eq 2 then begin
	catalog = a
	readcol, catalog, f='D, D', a, d
endif

;convert coordinate
if keyword_set(gal) then pathgal = 1b else pathgal = 0b	;path in gal or not?
fitsgal = strcmp(sxpar(hdr,'CTYPE1'), 'GL', 2, /fold_case)	;
if fitsgal eq pathgal then at=a & dt=d
if fitsgal and ~pathgal then glactc,a,d,2000,at,dt,1,/deg
if ~fitsgal and pathgal then glactc,at,dt,2000,a,d,2,/deg
adxy,hdr,at,dt,x,y	;convert wcs path to xy path on datacube

;convert polyline to spline
if keyword_set(spline) then spline_p,x,y,x,y	;if a spline resample is need

;resample path with interval of step
defstep=0.5	;unit pixel
if ~keyword_set(step) then step = defstep
if step le 0 then step = defstep
x = double(x) & y = double(y)
nodnum = min([n_elements(x),n_elements(y)])

px=x[0]
py=y[0]	;xy point on path
nod = 1	;xy of nod on the segment
path = [px, py]
while nod le nodnum-1 do begin
rest = sqrt((x[nod]-px)^2+(y[nod]-py)^2)
if rest ge step then begin ;the rest of current segment is enough for another step
	px += (x[nod]-px)/rest*step
	py += (y[nod]-py)/rest*step
	path = [path,px,py]
endif else begin	;not enough, find from the next segment or the next next one...except the last nod
	if nod eq nodnum-1 then begin	;the last nod, extend a little
		last = sqrt((x[nod]-px)^2+(y[nod]-py)^2)
		if last gt step/1e3 then begin
			px+=(x[nod]-px)/last*step
			py+=(y[nod]-py)/last*step
		path = [path, px, py]
		endif
		break
	endif
	while nod lt nodnum-1 do begin
		ptop2 = sqrt((x[nod+1]-px)^2 + (y[nod+1]-py)^2)
		if ptop2 lt step then nod++ else begin
			p1top2 = sqrt((x[nod+1]-x[nod])^2 + (y[nod+1]-y[nod])^2)
			dx21 = x[nod+1] - x[nod]
			dy21 = y[nod+1] - y[nod]
			crossx = x[nod]*y[nod+1]*dy21 - x[nod+1]*y[nod]*dy21 + dx21^2*px + dx21*dy21*py
			crossx /= dx21^2+dy21^2
			crossy = x[nod+1]*y[nod]*dx21 - x[nod]*y[nod+1]*dx21 + dx21*dy21*px + dy21^2*py
			crossy /= dx21^2+dy21^2
			crosstonext = sqrt((step^2 -(px-crossx)^2 -(py-crossy)^2) >0)
			px = crossx + dx21/p1top2*crosstonext
			py = crossy + dy21/p1top2*crosstonext
			path = [path,px,py]
			nod++
			break
		endelse
	endwhile
endelse
endwhile
path = reform(path,2,n_elements(path)/2)
pathx = (path[0,*])[*]
pathy = (path[1,*])[*]

;export the path
if keyword_set(track) then begin
	xyad,hdr,pathx,pathy,patha,pathd
	if fitsgal and ~pathgal then glactc,patha,pathd,2000,patha,pathd,2,/deg
	if ~fitsgal and pathgal then glactc,patha,pathd,2000,patha,pathd,1,/deg
	openw,lun,'slice.path',/get_lun
	printf,lun,transpose(string(patha)+' '+string(pathd))
	close,lun
	free_lun,lun
endif

;expand path to belt
stepnum = n_elements(pathx)
if keyword_set(width) then width=double(fix(width)) else width=1d
if width gt 1d then begin
	dy = shift(pathy, 1) - shift(pathy, -1)
	dy[0] = pathy[0]-pathy[1]
	dy[stepnum-1] = pathy[stepnum-2]-pathy[stepnum-1]
	dx = shift(pathx, 1) - shift(pathx, -1)
	dx[0] = pathx[0]-pathx[1]
	dx[stepnum-1] = pathx[stepnum-2]-pathx[stepnum-1]
	pa = atan(dy/dx)+(dx lt 0)*!pi+!pi/2
	dx = cos(pa)
	dy = sin(pa)
	pathx = pathx#replicate(1,width)+dx#(findgen(width)-(width-1)/2.)
	pathy = pathy#replicate(1,width)+dy#(findgen(width)-(width-1)/2.)
endif else width=1d

;resample the belt along the width
if keyword_set(pps) then begin
	pathx=congrid(pathx,n_elements(pathx[*,0]),2*n_elements(pathx[0,*]),/interp,/minus_one)
	pathy=congrid(pathy,n_elements(pathy[*,0]),2*n_elements(pathy[0,*]),/interp,/minus_one)
endif

;export the outline of the belt
if keyword_set(track) then begin
	if width gt 1d then begin
		xyad,hdr,pathx,pathy,patha,pathd
		if fitsgal and ~pathgal then glactc,patha,pathd,2000,patha,pathd,2,/deg
		if ~fitsgal and pathgal then glactc,patha,pathd,2000,patha,pathd,1,/deg
		ola = [patha[*,0], reverse(patha[*,width-1]), patha[0,0]]
		old = [pathd[*,0], reverse(pathd[*,width-1]), pathd[0,0]]
		openw,lun,'slice.belt',/get_lun
		printf,lun,transpose(string(ola)+' '+string(old))
		close,lun
		free_lun,lun
	endif
endif

;extract value from cube
if ~keyword_set(pps) then begin
	channum = sxpar(hdr, 'NAXIS3')
	slice = make_array(stepnum, width, channum, type=size(dat,/type))
	for i=0,channum-1 do slice[*,*,i] = interpolate(dat[*,*,i],pathx,pathy,missing=0)
	slice = mean(slice, dimension=2, /nan)
	slice = transpose(slice)
endif else begin
	slice = make_array(stepnum, width, type=size(dat,/type))
	slice = interpolate(dat,pathx,pathy,missing=0)
	slice = transpose(slice)
endelse

;output slice and resampled path
mkhdr,ppvhdr,slice

if ~keyword_set(pps) then begin
	sxaddpar,ppvhdr,'CTYPE1','VELOCITY'
	sxaddpar,ppvhdr,'CUNIT1','km/s'
	factor = abs(sxpar(hdr,'CDELT3') lt 100)?(1d):(1000d)
	sxaddpar,ppvhdr,'CRPIX1',sxpar(hdr,'CRPIX3')
	sxaddpar,ppvhdr,'CRVAL1',sxpar(hdr,'CRVAL3')/factor ;convert m/s to km/s
	sxaddpar,ppvhdr,'CDELT1',sxpar(hdr,'CDELT3')/factor
endif else begin
	sxaddpar,ppvhdr,'CTYPE1','POSITION'
	sxaddpar,ppvhdr,'CUNIT1','degree'
	sxaddpar,ppvhdr,'CRPIX1',1
	sxaddpar,ppvhdr,'CRVAL1',0d
	sxaddpar,ppvhdr,'CDELT1',abs(sxpar(hdr,'CDELT1'))/2
endelse

sxaddpar,ppvhdr,'CTYPE2','POSITION'
sxaddpar,ppvhdr,'CUNIT2','degree'
sxaddpar,ppvhdr,'CRPIX2',1
sxaddpar,ppvhdr,'CRVAL2',0d
sxaddpar,ppvhdr,'CDELT2',step*abs(sxpar(hdr,'CDELT1'))

if ~keyword_set(pps) then begin
	sxaddhist,'PV file: '+fitsfile,ppvhdr,/comment
	sxaddhist,'PV path:',ppvhdr,/comment
	for i=0,n_elements(a)-1 do sxaddhist,string(a[i])+' '+string(d[i]),ppvhdr,/comment
	widthnum = width*abs(sxpar(hdr,'CDELT1'))
	sxaddhist,'PV width: '+string(width)+' pixels'+string(widthnum*60)+' arcmin',ppvhdr,/comment

	fits_write,'pvslice.fits',slice,ppvhdr
endif else begin
	sxaddhist,'PP file: '+fitsfile,ppvhdr,/comment
	sxaddhist,'PP path:',ppvhdr,/comment
	for i=0,n_elements(a)-1 do sxaddhist,string(a[i])+' '+string(d[i]),ppvhdr,/comment
	widthnum = (width*2)*(abs(sxpar(hdr,'CDELT1'))/2)
        sxaddhist,'Pixel width: '+string(abs(sxpar(hdr,'CDELT1'))/2)+' arcdeg',ppvhdr,/comment
	sxaddhist,'PP width: '+string(width*2)+' pixels'+string(widthnum*60)+' arcmin',ppvhdr,/comment

	fits_write,'ppslice.fits',slice,ppvhdr
endelse

;print information
if ~keyword_set(pps) then print,'PV step: ',step,' * pixelsize' else print,'PP step: ',step,' * pixelsize'
if keyword_set(track) then print,'check *.path file for the info of the path'
if ~keyword_set(pps) then print,'PV width: ',width,' * pixelsize',widthnum*60,' arcmin' else begin
	print,'Pixel width: ',string(abs(sxpar(hdr,'CDELT1'))/2)+' arcdeg'
	print,'PP width: ',width*2,' * pixelsize',widthnum*60,' arcmin'
endelse
if keyword_set(track) then begin
  if width gt 1d then print,'check *.belt file for the info of the outline of the belt'
endif

end

