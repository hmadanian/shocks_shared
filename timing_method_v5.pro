; 
;
; PROCEDURE: timing_method_.pro 
; PURPOSE: To get bow shock  normal vector and velocity 
;          using four MMS spacecraft separation and dts 
;          in observation of the shock crossing
; KEYWORDS:
;         tplt: to plot the results (shifted and unshifted signals) on a tplot window
;         eps: save the tplot window as an .eps file
;
; NOTE: Results from this version are compared to the events in other papers
; Written by Hadi Madanian 
; 
; ## Accuracy notes ###  09/04/2023: output in the same ball park, but 
; ## sometime it matches other results sometimes it doesn't.
; ## also really depends on the time range for correlation.
; Testing events in other papers:
; Liu et al. 2020: FB boundary (https://iopscience.iop.org/article/10.3847/1538-4357/abb249)
; paper: 303   tm_v3: 325
; 
; Wang et al. 2020: flux rope speed (https://www.nature.com/articles/s41467-020-17803-3)
; paper: 119    tm_v3: 114 or 117
; 
; Yao et al. 2020: MHs (10.1029/2019JA026736)
; paper: 87 ± 6 km/s along [−0.72, −0.21, −0.65]
; tm_v3: 100.7 along   -0.76552780     0.080087267     -0.63839895
; tm_v5: 10/10/2024
;   Uploaded on github/shocks_shared
;
;


pro timing_method_v5, bdata = bdata, pdata = pdata, trange = trange, kvector = kvector, speed = speed, $
  window=window, data_rate = drate,  nobow = nobow, btype = btype, corrcoeff = corrcoeff, timedelays = timedelays, $
  eps = eps, tplt = tplt, matrix = matrix, tfinal = tfinal
  
@tplot_com

  current_window= !d.window
  if current_window ge 0 then begin
    currtpnm = tplot_vars.options.varnames
    currtptr  = tplot_vars.options.trange
  endif
  
  probes = ['1','2','3','4']
  if undefined(drate) then drate = 'brst'
  if undefined(bdata) then bdata = ['mms'+probes+'_fgm_b_gse_'+drate+'_l2_btot']
  if undefined(pdata) then pdata = ['mms'+probes+'_mec_r_gse']
  posdata = pdata
  if ~keyword_set(trange) then begin
    dprint,'Selecting time range for correlation from current tplot window: '
    if (current_window lt 0) then begin
      print, 'No active window found.... Error in multipoint calculations'
      stop
    endif
    print, 'Select time range for cross correlation"
    ctime, trange, npoints=2, /silent, /EXACT
    wait, 1
  endif
  tfinal = trange
  trt = [trange[0]-60, trange[0]+60]
 
  mustlod = 0
  mustlodpos = 0
  if ~tdexists(bdata[0], trange[0], trange[1]) then mustlod=1
  if ~tdexists(bdata[3], trange[0], trange[1]) then mustlod=1
  if ~tdexists(posdata[0], trt[0], trt[1]) then mustlodpos=1
  if ~tdexists(posdata[2], trt[0], trt[1]) then mustlodpos=1
  if (mustlod) then mms_load_fgm, trange=trt, probes = probes, data_rate= drate
  if (mustlodpos) then mms_load_state, trange = trt, probe = probes, datatypes = 'pos', /ephemeris_only
for ii=0,3 do tinterpol,posdata[ii], bdata[ii],/SPLINE,/NAN_EXTRAPOLATE 
posdata = posdata+'_interp'
store_data,'*intpd_shifted',/del

tb = trange
arrays = bdata
suffix = '_shifted'
dt_arr = dblarr(4) * 0.0
maxcorr = dt_arr
bn_mms1 = tsample(arrays[0], tb, times = tarray1, index = www)
dt = (tarray1[1]-tarray1[0]) 
bnt =n_elements(tarray1)
lags = lindgen(bnt) - ceil(bnt/2.0)
for ii = 1 , 3 do begin        ;;; Loop over each S/C pair
  curr_sc = arrays[ii]
  tinterpol, curr_sc , arrays[0], newname = curr_sc + '_intpd', /NEAREST_NEIGHBOR
  get_data, curr_sc + '_intpd', tarray2, dataarr
  result = dblarr(size(lags,/dim))
  for kk = 0 , n_elements(lags)-1 do begin
    shift_arr = shift(dataarr, lags[kk])
    bn_mms = shift_arr[www]
    result[kk] = CORRELATE(bn_mms1, bn_mms, /double)
  endfor
;  bn_mms = tsample(curr_sc + '_intpd', tb, times = tarray2)
;  result = C_CORRELATE(bn_mms1, bn_mms, lags, /double)
  maxcorr[ii] = max(result, maxind)
  dt_arr[ii] = -lags[maxind] * dt      ;;  f Finding the time lag
  get_data, curr_sc + '_intpd', data = dd
  store_data, curr_sc + '_intpd' + suffix, data = {x: (dd.x - dt_arr[ii]), y: dd.y}
  options, curr_sc + '_intpd' + suffix, 'labels', 'mms' + probes[ii]
endfor
; recreat MMS file variable for grouping later.
curr_sc = arrays[0]
copy_data, curr_sc, curr_sc + '_intpd' + suffix
copy_data, curr_sc, curr_sc + '_intpd'
options, curr_sc + '_intpd' + suffix, 'labels', 'mms' + probes[0]
options, curr_sc + '_intpd', 'labels', 'mms' + probes[0]

possc = dblarr(3,4) * 0.0
tcros1 = dblarr(4) * 0.0
t1 = tarray1[ceil(bnt/2)]
tcros1 = dt_arr + t1
Tmat = dt_arr[1:-1] ; [tcros1[1] - t1, tcros1[2] - t1, tcros1[3] - t1]
for ii=0,n_elements(tcros1)-1 do possc[*,ii] = data_cut(posdata[ii],tcros1[ii]) ;mms3_fgm_r_gse_brst_l2_vec, or mms1_mec_r_gse
inds = where(possc eq 0.0, cnt)
if cnt gt 0 then begin
  print, 'Position vector is empty: Somethin not right'
  stop
endif
Dmat = [[possc[*,1]-possc[*,0]], [possc[*,2]-possc[*,0]], [possc[*,3]-possc[*,0]]]
m =LA_LINEAR_EQUATION(Dmat, Tmat)
kvector = m /norm(m, /double,LNORM=2)
Vwave = m /norm(m, /double,LNORM=2) / norm(m, /double,LNORM=2)   ; Vector
speed = 1/norm(m, /double,LNORM=2)
corrcoeff = maxcorr[1:-1]
timedelays = dt_arr[1:-1]
inds = where(maxcorr[1:-1] lt 0.85,cnt)
if cnt gt 0 then begin
  print, '!!!!=======---- Correlation coefficient too low < 0.85, bad correlation, maybe choose a differnet period'
  print, 'between MMS1 and', probes[inds+1]
endif

dsize = get_screen_size()
if keyword_set(tplt) then begin
  if keyword_set(window) then Dwin = window else begin
    window, /free, xsize=dsize[0]/2., ysize=dsize[1]*2.5/3.,xpos=dsize[0]/4., ypos=dsize[1]/3.
    Dwin = !d.window
  endelse
  wset,Dwin
  stornm = 'unshifted_original'
;  store_data, stornm, /del
  store_data, stornm, data = bdata
  options, stornm ,'colors', [0, 210, 1, 135]
  options, stornm, 'labels', ['mms'+probes]
  options, stornm, 'labflag', -1
  options, bdata,'ystyle',1
  options, bdata,'yrange',[0,0]

  stornm = 'shifted_signals'
;  store_data, stornm, /del
  newarrays = bdata + '_intpd' + suffix
  store_data, stornm, data = [ newarrays ], dlim={constant:0} ;['mms2_Bvec_brst_filt_ntt_x_shift','mms3_Bvec_brst_filt_ntt_x_shift','mms4_Bvec_brst_filt_ntt_x_shift','mms1_Bvec_brst_filt_ntt_x']
  options, stornm ,'colors', [0, 210, 1, 135]
  options, stornm, 'labels', ['mms'+probes]
  options, stornm, 'labflag', -1
  options, newarrays,'ystyle',1
  options, newarrays,'yrange',[0,0]
  options,'*_intpd_shifted','ystyle',1
  tplot,['unshifted_original', 'shifted_signals'], window=Dwin
  timebar, tb, linestyle=2
endif
;wait, 0.2
store_data,posdata,/del

print, 'Correlation Coefs: ', corrcoeff
print, 'Time delays (s): ', timedelays
print, 'k_vec in gse: ', kvector
print, 'speed: ', speed
print, 'final time range:', time_string(tfinal)
if keyword_set(matrix) then begin
  print, 'k_vec in NCB: ', matrix ## kvector
endif

if keyword_set(eps) then begin
  popen, 'timing_method_figure', /encapsulated, xsize=8, ysize=4, unit='inches'
  tplot,['unshifted_original', 'shifted_signals']
  timebar, tb, linestyle=2
  xyouts, 0.3 , 0.4, 'Correlation Coefs: ' + strtrim(corrcoeff, 2), /NORMAL, color=0
  xyouts, 0.3 , 0.5, 'Time delays (s): ' + strtrim(timedelays, 2), /NORMAL, color=0
  xyouts, 0.3 , 0.6, 'k_vec in gse: ' + strtrim(kvector, 2), /NORMAL, color=0
  xyouts, 0.3 , 0.7, 'speed: ' + strtrim(speed, 2) + ' km/s', /NORMAL, color=0
  pclose
endif 


;store_data, '*_intpd*',/del

if (current_window ge 0) then tplot, currtpnm, window=current_window, trange = currtptr































end