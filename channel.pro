;Purpose: dervie a new fitsfile with desired channels from the original fitsfile
;Usage: channel, fitsfile, v1, v2, num, rmsfile=rmsfile
;Input:
;  fitsfile: the original fitsfile WITHOUT suffix
;  v1: the CENTER velocity of the first desired channprintel
;  v2: the CENTER velocity of the last desired channel
;  num: is the number of the desired channels
;Input keyword:
;  rmsfile: the rms fitsfile used to set the value under 3sigma to 0 in each channel
;Caution: 
;  be aware that the integrated range is set as 1 km/s, 
;    which is from (CENTER velocity - 0.5 km/s) to (CENTER velocity + 0.5 km/s)
;  the unit of v1 and v2 is km/s
;  the rms fitfile is needed to set the keyword rmsfile
;History:
;  Completed. 22 Nov 2015.
;  Update the Head Info of resulting fitsfile. 7 Dec 2015. 
;Copyright: modified by fxong@PMO

pro v2c,hdr,v,c
        if n_params() lt 3 then begin
                print, 'Syntax - V2C, hdr, v, c'
                return
        endif
        nc = sxpar(hdr,'NAXIS3')
        v0 = sxpar(hdr,'CRVAL3')
        c0 = sxpar(hdr,'CRPIX3')
        dv = sxpar(hdr,'CDELT3')
        c = (v-v0)/dv+c0-1
end

;hrebinv
;by ShaoboZhang
;History:
;Jul,23,2015,v1.0

pro hrebinv, oldcub, oldhdr, newcub, newhdr, velo1, velo2, nchannel
;Rebin velocity axis for channel map
;the given velo1 and velo2 is the CENTER velocity of the first and last channel
;Input:
;	oldcub, oldhdr: the data and header of old fits
;	velo1, velo2: scalars indicate the velocity range
;	nchannel: number of channel of result cube
;Output:
;	newcub, newhdr: the data and header for new fits
;Usage: hrebinv, oldcub, oldhdr, newcub, newhdr, v1, v2, nc
	if n_params() lt 2 then begin
	    print, 'Syntax - HREBINV, oldcub, oldhdr, newcub, newhdr, velocity1, velocity2, num_channels'
	    return
	endif
	v2c,oldhdr,velo1,c1
	v2c,oldhdr,velo2,c2
	width = float(abs(c1-c2))/(nchannel-1)
	c = min([c1,c2])-width/2 + findgen(nchannel+1)*width
	olddim = size(oldcub,/dimension)
	newcub = fltarr(olddim[0],olddim[1],nchannel)

	for i=0,nchannel-1 do begin
		if c[i] lt -0.5 or c[i+1] gt olddim[2]-0.5 then begin
			newcub[*,*,i] = !values.f_nan
			print, 'Warning - given velocity range exceed the old datacube'
			continue
		end
		intc1 = ceil(c[i]+0.5)
		fltc1 = intc1-c[i]-0.5
		intc2 = floor(c[i+1]-0.5)
		fltc2 = c[i+1]-intc2-0.5

		if fltc1 eq 0 then head = 0 else head = (oldcub[*,*,intc1-1] * fltc1)
		if fltc2 eq 0 then rear = 0 else rear = (oldcub[*,*,intc2+1] * fltc2)
		case intc1-intc2 of
			0:begin
				newcub[*,*,i] = head + oldcub[*,*,intc1] + rear
			end
			1:begin
				newcub[*,*,i] = head + rear
			end
			2:begin
				newcub[*,*,i] = head /fltc1 *(fltc1+fltc2-1)
			end
			else:begin
				newcub[*,*,i] = head + total(oldcub[*,*,intc1:intc2],3) + rear
			end
		endcase
	endfor
	;newcub *= sxpar(oldhdr, 'CDELT3')	;convert sum(Tmb) to integrated intensity
	newcub *= sxpar(oldhdr, 'CDELT3')/1000d
	newhdr = oldhdr
	sxaddpar, newhdr, 'NAXIS3', nchannel
	sxaddpar, newhdr, 'CRPIX3', 1
	sxaddpar, newhdr, 'CRVAL3', velo1
	sxaddpar, newhdr, 'CDELT3', (velo2-velo1)/(nchannel-1)
end

;Main Procedure
pro channel, fitsfile, v1, v2, num, rmsfile=rmsfile

if n_params() lt 3 then begin
  print,'Syntax - CHANNEL, fitsfile, velocity1, velocity2, num_channels, [rmsfile= ]'
  return
endif
if ~file_test(fitsfile+'.fits') then begin
  print,'Error: Fits file does not exist! Input fits file without suffix.'
  return
endif
if keyword_set(rmsfile) then begin
  if ~file_test(rmsfile+'.fits') then begin
    print,'Error: Rms file does not exist! Input rms file without suffix.'
    return
  endif
endif

fits_read, fitsfile+'.fits', dat, hdr
factor=abs(sxpar(hdr,'CDELT3') lt 100)?(1d):(1000d)
hrebinv, dat, hdr, newdat, newhdr, v1*factor, v2*factor, num
sxaddpar,newhdr,'CRVAL3',sxpar(newhdr,'CRVAL3')/factor
sxaddpar,newhdr,'CDELT3',sxpar(newhdr,'CDELT3')/factor
sxdelpar,newhdr,['NAXIS','CTYPE','CRVAL','CDELT','CRPIX','CROTA']+'4' ;MWISP fits
sxdelpar,newhdr,['ORIGIN','DATE','RMS','HISTORY']                     ;MWISP fits
sxaddhist,'Number of channels: '+string(sxpar(newhdr,'NAXIS3')),newhdr,/comment
sxaddhist,'Central Velocity:',newhdr,/comment
for i=0, sxpar(newhdr,'NAXIS3')-1 do begin
  sxaddhist,string(sxpar(newhdr,'CRVAL3')+sxpar(newhdr,'CDELT3')*i),newhdr,/comment
endfor
sxaddhist,'Velocity in km/s',newhdr,/comment
sxaddhist,'Integrated Range: (Central V - 0.5 km/s) to (Central V + 0.5 km/s)',newhdr,/comment
sxaddhist,'Channel file: '+fitsfile+'.fits',newhdr,/comment

if keyword_set(rmsfile) then begin
  fits_read,rmsfile+'.fits',rdat,rhdr
  dev=sxpar(hdr,'CDELT3')/factor
  Irms=rdat*sqrt(dev)
  for k=0, sxpar(newhdr,'NAXIS3')-1 do begin
    for i=0, sxpar(newhdr,'NAXIS1')-1 do begin
      for j=0, sxpar(newhdr,'NAXIS2')-1 do begin
        if newdat[i,j,k] lt 3*Irms[i,j] then begin
          newdat[i,j,k]=0
        endif
      endfor
    endfor
  endfor
  fits_write, fitsfile+'_cha_rm.fits', newdat, newhdr
endif else begin
  fits_write, fitsfile+'_cha.fits', newdat, newhdr
endelse

print,'the Number of the desired channels:', sxpar(newhdr,'NAXIS3')
print,'the Central Velocity of each channel:'
for i=0, sxpar(newhdr,'NAXIS3')-1 do begin
  print,sxpar(newhdr,'CRVAL3')+sxpar(newhdr,'CDELT3')*i,' km/s'
endfor
print,'the Interval of the desired channels:', sxpar(newhdr,'CDELT3'),' km/s'
print,'the Integrated Range is from (Central V - 0.5 km/s) to (Central V + 0.5 km/s).'

end
