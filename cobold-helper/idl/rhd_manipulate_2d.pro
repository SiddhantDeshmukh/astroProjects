; pro rhd_manipulate
;
;------------------------------------------------------------------------------
;+
; NAME:
;   rhd_manipulate_2d ('radiation_hydrodynamics_manipulate')
;
; PURPOSE:
;
; CATEGORY:
;
; CALLING SEQUENCE:
;   @rhd_manipulate
;
; INPUTS:
;
; OUTPUTS:
;
; KEYWORDS:
;
; COMMON BLOCKS:
;
; ROUTINES:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;   Before first use type
;     delvar, eos
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   17. 5.99 Written by Bernd Freytag, Astronomical Observatory Copenhagen
;    8. 6.00 (B.F.) Add points to model
;   30. 1.01 (B.F.)
;   11.12.01 (B.F.) Subroutine introduced
;   04.02.21 (S.A.D.) Changed 'name_quc' and 'name_species' for CHEM module
;   05.02.21 (S.A.D.) New routine to shrink 3D model to 2D by eliminating 'y'
;------------------------------------------------------------------------------

; Solar model
model='d3t57g44mm20'& modelident='n03'  ; relaxed 3D solar model
modelhome = '/home/sd/cobold/runs/start/'

outmodelhome = '/home/sd/cobold/runs/start/'
outmodel    =modelhome + model & outmodelident='' + ''
outmodel = 'd2gt57g44mm20'
outmodelident = ''

; Metal-poor model
;model='d3t57g44mm20n03' & modelident=''  ; 3D metal-poor model
;modelhome = '/home/sd/cobold/runs/metal-poor/start/'
;outmodelhome = '/home/sd/cobold/runs/metal-poor/start/'

;outmodel = 'd2t57g44mm20n03'
;outmodelident = '.00'
; No disk, normal model with CHEM quc
;parfile = modelhome + model + modelident + '.par'
modelfile = modelhome + model + modelident + '.start' & ndata=0
outmodelfile = outmodelhome + outmodel + outmodelident + '.start'
; --- Read requested model number ---
print, 'Read initial model from ' + modelfile
ful=uio_dataset_rd(modelfile, n=ndata, /HeadRead)
;eosbox, ful, eos=EOS, opa=1, out_q=0, ierror=ierror & help,ierror
;
; --- Physical units ---
unit_time='s'
unit_x1  ='cm'
unit_x2  =unit_x1
unit_x3  =unit_x1
unit_rho ='g/cm^3'
unit_v   ='cm/s'
unit_e   ='erg/g'
unit_B   ='G*sqrt(4pi)'
;
; --- Time information ---
itime=ful.z.itime*0L           ; --- *0L
time =ful.z.time*0.0           ; --- *0.0
dtime=ful.dtime                ; --- *0.5 *0.7 *0.95
time_out_full_last=-1.0E9
time_out_mean_last=time_out_full_last
;
; --- Put data from input structure into individual arrays ---
;     indexstart='central'
print,'xb2:', ful.z.xb2[0:1]
rhd_manipulate_parsebox, ful, indexstart='box', $
      m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, Bb1=Bb1,Bb2=Bb2,Bb3=Bb3
delvar,quc
print,'xb2:', xb2[0:1]
if (n_elements(quc) ne 0) then nquc=n_elements(quc[0,0,0,*])
head=ful.head
;delvar,ful
;
; Interpolate onto 2D
; Here 'num2' doesn't matter so long as it's not zero or n2-m2+1 since
; the procedure just copies the 'xc1' bit and removes all but the first
; part (slice at '0')
rhd_manipulate_interpolate_2d, 0, 1, 0, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3

help, itime, time, dtime, v1
help, m1,n1, m2,n2, m3,n3
help, xc2,xb3
print, 'minmax(rho,ei,v1,v2,v3)'
print, minmax(rho) & print, minmax(ei) & print, minmax(v1) & print, minmax(v2) & print, minmax(v3)
rhomin=min(rho)
;
; Initialise CHEM input
;nquc=8
; Better to make 'name_species' and then add 'Number density of' as a prepend string after the call to 'co_initialabu()'
; This is just written out explicitly for now
;name_quc='Number density of '+['H','H2','C','O','CO','CH','OH','metal']
;name_species=['H','H2','C','O','CO','CH','OH','metal']  ; Same species as 'name_quc', passed to init abu stuff
;unit_quc=replicate('cm^-3',nquc)
;quc=fltarr(szq[0],szq[1],szq[2],8)
;quc=fltarr(n1-m1+1,n2-m2+1,n3-m3+1,8)

; metal abundance (as fraction of H density)
;metal=0.0851

;atomdata=co_loadatomdata()

; Calculate initial number density for each grid cell
;for i3=0,n3-m3 do begin&$
;  for i2=0,0 do begin&$
;    for i1=0,n1-m1 do begin&$
;    quc(i1,i2,i3,*)=float(co_initialabu(rho(i1,i2,i3),atomdata,name_species,metal=metal))&$
;    endfor&$
;  endfor&$
;endfor&$

;
; --- End processing ----------------------------------------------------------
;
rhd_manipulate_assure3d, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3
;
help, m1,n1, m2,n2, m3,n3
print, 'minmax(rho,ei,v1,v2,v3)'
print, minmax(rho) & print, minmax(ei) & print, minmax(v1) & print, minmax(v2) & print, minmax(v3)
;
; --- Model information ---
fileform='unformatted' & fileconv='ieee_4' & double_flag=0
version    ='rhd_manipulate 2017-10-24'
description=[head.description, $
             'Based on: ' + modelfile, $
             'Add points on all sides: 40: 317^3->397^3 points, rhofac=0.95']

if (n_elements(head.history) gt 13) then begin &$
  history    =[head.history[0:8], $
               '...', $
               head.history[12:*], $
               'Processed by rhd_manipulate: ' + strtrim(systime(0),2)] &$
endif else begin &$
  history    =[head.history, $
               'Processed by rhd_manipulate: ' + strtrim(systime(0),2)] &$
endelse
;
; stop ; ----------------------------
;
print, '> Make box'
rhd_manipulate_makebox, ful1, $
     m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, Bb1=Bb1,Bb2=Bb2,Bb3=Bb3
;
; stop ; ----------------------------
;
; --- Write data into file ---
; rhd_wrboxdata, 'open,model,box,endmo,close', $
;    -1, ncout, file=outmodelfile, form=fileform, conv=fileconv, double_flag=double_flag, $
;    ndes=n_elements(description), description=description, $
;    nhis=n_elements(history), history=history, $
;    version=version, $
;    unit_time=unit_time, $
;    unit_x1=unit_x1, unit_x2=unit_x2, unit_x3=unit_x3, $
;    unit_rho=unit_rho, unit_v=unit_v, unit_e=unit_e, $
;    unit_quc=unit_quc, name_quc=name_quc, unit_B=unit_B, $
;    m1=m1, m2=m2, m3=m3, n1=n1, n2=n2, n3=n3, $
;    itime=itime, time=time, dtime=dtime, $
;    out_time_full_last=time_out_full_last, out_time_mean_last=time_out_mean_last, $
;    xc1=xc1, xc2=xc2, xc3=xc3, xb1=xb1, xb2=xb2, xb3=xb3, $
;    rho=rho, v1=v1, v2=v2, v3=v3, ei=ei, $
;    quc=quc, Bb1=Bb1, Bb2=Bb2, Bb3=Bb3
;
rhd_wrboxdata, 'open,model,box,endmo,close', $
   -1, ncout, file=outmodelfile, form=fileform, conv=fileconv, $
   ndes=n_elements(description), description=description, $
   nhis=n_elements(history), history=history, $
   version=version, $
   unit_time=unit_time, $
   unit_x1=unit_x1, unit_x2=unit_x2, unit_x3=unit_x3, $
   unit_rho=unit_rho, unit_v=unit_v, unit_e=unit_e, $
   unit_quc=unit_quc, name_quc=name_quc, $
   unit_B=unit_B, $
   m1=m1, m2=m2, m3=m3, n1=n1, n2=n2, n3=n3, $
   itime=itime, time=time, dtime=dtime, $
   out_time_full_last=time_out_full_last, out_time_mean_last=time_out_mean_last, $
   xc1=xc1, xc2=xc2, xc3=xc3, xb1=xb1, xb2=xb2, xb3=xb3, $
   rho=rho, v1=v1, v2=v2, v3=v3, ei=ei, $
   Bb1=Bb1, Bb2=Bb2, Bb3=Bb3, quc=quc

;sta=uio_dataset_rd(outmodelfile)
sta=uio_dataset_rd(outmodelfile, /HeadRead)
;
; delvar, rho,ei,v1,v2,v3,quc
; print,([1.0E-04, 5.0E-04, 2.0E-03]*4.0)^(-0.5) * (!con.Rsun/!con.au)
; print,([2.0E-04, 6.0E-04,18.0E-04]*4.0)^(-0.5) * (!con.Rsun/!con.au)
;
;end ; rhd_manipulate
