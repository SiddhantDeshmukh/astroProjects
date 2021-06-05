; pro rhd_manipulate
;
;------------------------------------------------------------------------------
;+
; NAME:
;   rhd_manipulate ('radiation_hydrodynamics_manipulate')
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
;   04.02.21 (S.A.D) Changed 'name_quc' and 'name_species' for CHEM module
;------------------------------------------------------------------------------

; Change initial model
model='d2t59g40mm20' & modelident='.171'
modelhome = '/home/sd/cobold/runs/full/d2t59g40mm20/'
outmodelhome = '/home/sd/cobold/runs/start/'
;
outmodel    = model & outmodelident='' + 'c269co072cu12'

; No disk, normal model with CHEM quc
parfile = modelhome + model + modelident + '.par'
modelfile = modelhome + model + modelident + '.full' & ndata=0
outmodelfile = outmodelhome + outmodel + outmodelident + '.start'
;
; ; --- Read EOS file ---
; par=uio_struct_rd(parfile)
; if (n_elements(EOS) eq 0) then begin &$
;   tabinter_rdcoeff, eosdisk + par.eosfile, EOS &$
;   dfopta, opadisk + par.opafile &$
; endif
; f0=!con.sigmasb*par.Teff^4
; ;
; ; --- Dust parameters ----
; ms =140.71*!con.muAtom &$                         ; (2.0*24.31+1.0*28.09+4.0*16.0)
; dustfrac=6.0324E-04/(2.0*24.31)*ms/!con.muAtom &$ ; AGS abu, Mg2SiO4&$
; ;
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
;name_quc =['rho_monomer', 'rho_dustr1']
;unit_quc =['g/cm^3'     , 'g/cm^3'    ]
;name_quc =['n_forsterite', 'n_corundum']
;unit_quc =['1/cm^3'      , '1/cm^3'      ]
;name_quc =['rho_monomer', 'rho_dustr1', 'rho_dustr2', 'rho_dustr3', 'rho_dustr4', 'rho_dustr5']
;unit_quc =['g/cm^3'     , 'g/cm^3'    , 'g/cm^3'    , 'g/cm^3'    , 'g/cm^3'    , 'g/cm^3'    ]
;name_quc =['rho_monomer', 'rho_dustr1', 'rho_dustr2', 'rho_dustr3', 'rho_dustr4', 'rho_dustr5', 'rho_dustr6', 'rho_dustr7']
;unit_quc =['g/cm^3'     , 'g/cm^3'    , 'g/cm^3'    , 'g/cm^3'    , 'g/cm^3'    , 'g/cm^3'    , 'g/cm^3'    , 'g/cm^3'    ]
;name_quc =['']
;unit_quc =['']
;unit_B   ='G'
unit_B   ='G*sqrt(4pi)'
;
; --- Time information ---
itime=ful.z.itime*0L           ; --- *0L
time =ful.z.time*0.0           ; --- *0.0
dtime=ful.dtime                ; --- *0.5 *0.7 *0.95
;dtime=1.0E+02
time_out_full_last=-1.0E9
;time_out_full_last=ful.time_out_full_last
;time_out_mean_last=ful.time_out_mean_last
;time_out_mean_last=ful.time_out_full_last
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
help, itime, time, dtime, v1
help, m1,n1, m2,n2, m3,n3
print, 'minmax(rho,ei,v1,v2,v3)'
print, minmax(rho) & print, minmax(ei) & print, minmax(v1) & print, minmax(v2) & print, minmax(v3)
rhomin=min(rho)
;
; --- Start processing --------------------------------------------------------

; --- Damp velocities outside certain radius ---
;radius=combox('xc3a',ful)
;radius=radius-min(radius)
;radius=radius/max(radius)
;dampradius=sqrt(exp(-(radius/0.80)^16))
;radius=combox('radiusnorm',ful)
;dampradius=(radius/95.0E+12)^4/(radius/95.0E+12 + 1.0)^4
;dampradius=(1.0-cosinebell(radius, 0.25, 0.05))
;dampradius=sqrt(exp(-(radius/0.50)^16))
;dampradius=sqrt(exp(-(radius/0.43)^16))
;v1=v1*dampradius
;v2=v2*dampradius
;v3=v3*dampradius

; --- Remove all horizontal fluctuations (in base quantities) ---
;v1=v1*0.0
;v2=v2*0.0
;v3=v3*0.0
;rhoavg  =avg(avg(rho   ,0),0)
;rhoeiavg=avg(avg(rho*ei,0),0)
;rho=rebin(reform(rhoavg         , 1,1,n3-m3+1), n1-m1+1, n2-m2+1, n3-m3+1)
;ei =rebin(reform(rhoeiavg/rhoavg, 1,1,n3-m3+1), n1-m1+1, n2-m2+1, n3-m3+1)
;qucavg  =avg(avg(quc[*,*,*,0],0),0)
;quc[*,*,*,0]=rebin(reform(qucavg, 1,1,n3-m3+1), n1-m1+1, n2-m2+1, n3-m3+1)
;qucavg  =avg(avg(quc[*,*,*,1],0),0)
;quc[*,*,*,1]=rebin(reform(qucavg, 1,1,n3-m3+1), n1-m1+1, n2-m2+1, n3-m3+1)

; --- Damp fluctuations above/below certain depth (given by density) ---
;vmask=((10.0^(-7.4)/rho) < 1.0)^3
;vmask=((10.0^(-4.7)/rho) < 1.0)^4
;vmask=cosinebell((rho-max(rho))/(max(rho)-min(rho))*0.5, 0.2)
;vmask=((rho/10.0^(-10.5)) < 1.0)^3
;vmask=((rho/10.0^(-9.0)) < 1.0)^3
;vmask=((rho/10.0^(-8.0)) < 1.0)^3
;vmask=((rho/10.0^(-4.0)) < 1.0)^3
;vmask=((rho/10.0^(-4.4)) < 1.0)^3
;v1=v1*vmask
;v2=v2*vmask
;v3=v3*vmask
;rho0    =rho
;rhoeiavg=rebin(reform(avg(avg(rho*ei,0),0),1,1,n3),n1,n2,n3)
;rhoavg  =rebin(reform(avg(avg(rho   ,0),0),1,1,n3),n1,n2,n3)
;rho     = rhoavg  +vmask*(rho    -rhoavg  )
;ei      =(rhoeiavg+vmask*(rho0*ei-rhoeiavg))/rho
;qucavg  =rebin(reform(avg(avg(quc,0),0),1,1,n3,nquc),n1,n2,n3,nquc)
;quc[*,*,*,0]=rho/rhoavg*qucavg[*,*,*,0]+vmask*(quc[*,*,*,0]-rho/rhoavg*qucavg[*,*,*,0])
;quc[*,*,*,1]=rho/rhoavg*qucavg[*,*,*,1]+vmask*(quc[*,*,*,1]-rho/rhoavg*qucavg[*,*,*,1])

; --- Remove points on all sides ---
;addn=42 & rhd_manipulate_rempoints, addn, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3
;addn=12 & rhd_manipulate_rempoints, addn, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3

; --- Shift axes ---
;xc1=xc1-xb1[0]
;xb1=xb1-xb1[0]
;xc2=xc2-xb2[0]
;xb2=xb2-xb2[0]
;xc3=xc3-xb3[0]
;xb3=xb3-xb3[0]

; --- Interpolate onto other grid ---
;rhd_manipulate_interpolate, 800, 800, 400, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='1'
;rhd_manipulate_interpolate, 560, 560, 280, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='1'
;rhd_manipulate_interpolate, 400, 400, 200, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='1'
;rhd_manipulate_interpolate, 350, 350, 140, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='1'
;rhd_manipulate_interpolate, 700, 700, 140, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='1'
;rhd_manipulate_interpolate, 1400, 1400, 140, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='1'
;rhd_manipulate_interpolate, 197, 197, 197, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='central'
;rhd_manipulate_interpolate, 281, 281, 281, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='central'
;rhd_manipulate_interpolate, 317, 317, 317, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='central'
;rhd_manipulate_interpolate, 361, 361, 361, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='central'
;rhd_manipulate_interpolate, 401, 401, 401, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='central'
;rhd_manipulate_interpolate, 425, 425, 425, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='central'
;rhd_manipulate_interpolate, 479, 479, 479, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='central'
;rhd_manipulate_interpolate, 2*(n1-m1+1), 0, 0, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='box'
;rhd_manipulate_interpolate, 300, 300, 0, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, indexstart='box'
;rhd_manipulate_interpolate, 224, 224, 0, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, indexstart='box'
;rhd_manipulate_interpolate, 448, 448, 0, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, indexstart='box'
;print,'xb2:', xb2[0:1]
;rhd_manipulate_interpolate, 0, (n1-m1+1), 0, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, indexstart='box'
;rhd_manipulate_interpolate, 0, 180, 0, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='box'
;dx3=(max(xb1)-min(xb1))/((n1-m1+1)*2.0)
;n3new=fix((max(xb3)-min(xb3))/dx3)+2 & print, 'n3 -> n3new', n3, n3new
;xb3new=findgen(n3new)*dx3+min(xb3)
;rhd_manipulate_interpolate, 0, 0, xb3new, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, indexstart='box', x3_command='newxb3'
;xb3new=findgen(163)*14.0E+05
;rhd_manipulate_interpolate, 0, 0, xb3new, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, indexstart='1', x3_command='newxb3'
;rhd_manipulate_interpolate, 1000, 0, 400, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, indexstart='box'

;xfac=10.0^0.4
;xfac=1.422222
;xfac=2.0
;xfac=2.5
;xfac=0.8
;xfac=0.7
;xfac=900.0/1020.0
;xfac=630.0/570.0
;xfac=580.0/540.0
;xfac=360.0/352.0
;xc1=xc1*xfac & xb1=xb1*xfac
;xc2=xc2*xfac & xb2=xb2*xfac
;xc3=xc3*xfac & xb3=xb3*xfac
;rho=rho/xfac
;v3=v3/xfac^2
;v3=v3/xfac
;v1=v1*xfac
;v2=v2*xfac
;xc2=xc2*300.0 & xb2=xb2*300.0
;xc2=xc2*0.5 & xb2=xb2*0.5
;print,'xb2:', xb2[0:1]
;xc2=xc1 & xb2=xb1
;print,'xb2:', xb2[0:1]

;xb3(1:n3-m3)=0.5 * (xc3(0:n3-m3-1) + xc3(1:n3-m3  ))
;xb3(n3-m3+1)=xb3(n3-m3) + (xb3(n3-m3)-xb3(n3-m3-1))^2/(xb3(n3-m3-1)-xb3(n3-m3-2))

; --- Add points ---
; addn=40 & rhd_manipulate_addpoints, addn, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, rhofac=0.95

;n3shift=m3-1
;m3=m3-n3shift
;n3=n3-n3shift

;dx1=min(xb1) & xc1=xc1-dx1 & xb1=xb1-dx1
;dx2=min(xb2) & xc2=xc2-dx2 & xb2=xb2-dx2
;dx3=max(xb3) & xc3=xc3-dx3 & xb3=xb3-dx3
;dx3=min(xb3) & xc3=xc3-dx3 & xb3=xb3-dx3

; --- Remove points (one sided) ---
;subn=18
;n3new=n3-subn
;xc3=xc3[    0:n3new-m3]
;xb3=xb3[    0:n3new-m3+1]
;v1 =v1[ *,*,0:n3new-m3]
;v2 =v2[ *,*,0:n3new-m3]
;v3 =v3[ *,*,0:n3new-m3]
;ei =ei[ *,*,0:n3new-m3]
;rho=rho[*,*,0:n3new-m3]
;quc=quc[*,*,0:n3new-m3,*]
;n3=n3new


; --- Replicate the box horizontally ---
;rhd_manipulate_repeat, 2, 2, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc

; --- Add global rotation to model velocities ---
;rhd_manipulate_rotate, axis, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3
;
; --- Add random velocities ---
;v3randomfac=1.0E-01 & v3=v3*(1.0+v3randomfac*randomn(seed, n1-m1+1, n2-m2+1, n3-m3+1))
;v1randomfac=1.0E+01 & v1=v1+v1randomfac*randomn(seed, n1-m1+1, n2-m2+1, n3-m3+1)
;v2randomfac=1.0E+02 & v2=v2+v2randomfac*randomn(seed, n1-m1+1, n2-m2+1, n3-m3+1)
;v3randomfac=1.0E+02 & v3=v3+v3randomfac*randomn(seed, n1-m1+1, n2-m2+1, n3-m3+1)

; --- Remove average velocity ---
;rhoavg=avg(avg(rho,0),0)
;rhov1avg=avg(avg(rho*v1,0),0) & v1avg=rhov1avg/rhoavg & v1avg=rebin(reform(v1avg, 1,1,n3-m3+1), n1-m1+1, n2-m2+1, n3-m3+1)
;v1=v1-v1avg
;rhov2avg=avg(avg(rho*v2,0),0) & v2avg=rhov2avg/rhoavg & v2avg=rebin(reform(v2avg, 1,1,n3-m3+1), n1-m1+1, n2-m2+1, n3-m3+1)
;v2=v2-v2avg
;rhoavg=avg(avg(rho,0),0) & rhov3avg=avg(avg(rho*v3,0),0) & v3avg=rhov3avg/rhoavg & v3=v3-v3avg & v3avg=rebin(reform(v3avg, 1,1,n3-m3+1), n1-m1+1, n2-m2+1, n3-m3+1)
;v3=v3-v3avg
;
;vfac=0.1 & v1=vfac*v1 & v2=vfac*v2 & v3=vfac*v3

; --- Create dust array ---
;quc=replicate(0.0, n1-m1+1, n2-m2+1, n3-m3+1, 2)
;quc[*,*,*,0]=rho*dustfrac
;quc[*,*,*,1]=rho*dustfrac*1.0E-10

; --- Remove dust above certain height, reduce amount of momomers ---
;quc[*,*,90:*,0  ]=quc[*,*,90:*,0]*1.0E-03
;quc[*,*,90:*,1:*]=0.0
;quc[*,*,*,0]=0.0
;quc[*,*,130:*,*]=1.0
;quc[*,*,0:180,1]=quc[*,*,0:180,1] < (0.999*dustfrac*rho[*,*,0:180])
;quc[*,*,0:180,0]=((dustfrac*rho[*,*,0:180]) - quc[*,*,0:180,1])
;r=       rebin(reform([xc1^2], n1-m1+1,       1,       1), n1-m1+1, n2-m2+1, n3-m3+1)
;r=     r+rebin(reform([xc2^2],       1, n2-m2+1,       1), n1-m1+1, n2-m2+1, n3-m3+1)
;r=sqrt(r+rebin(reform([xc3^2],       1,       1, n3-m3+1), n1-m1+1, n2-m2+1, n3-m3+1))
;ind=where(r gt max(xb3))
;r=0.0
;quc[ind]=(1.0E-10*dustfrac)*rho[ind]
;quc[ind+((n1-m1+1L)*(n2-m2+1L)*(n3-m3+1L))]=(1.0E-10*dustfrac)*rho[ind]
;ind=0
;quc[*,*,*,0]=(rho*dustfrac) < quc[*,*,*,0]
;quc[*,*,*,1]=0.0 > ((rho*dustfrac-quc[*,*,*,0]) < quc[*,*,*,1])

; --- Distribute dust over more bins ---
;qucold=quc
;nqucold=nquc
;nquc=8
;quc=fltarr(n1-m1+1, n2-m2+1, n3-m3+1, nquc)
;; ... Copy monomers ...
;quc[*,*,*,0]=qucold[*,*,*,0]
;; ... Distribute dust equally over first nqf bins ...
;nqf=5 & quc[*,*,*,1:nqf]=rebin(reform(qucold[*,*,*,1]/float(nqf), n1-m1+1, n2-m2+1, n3-m3+1, 1), n1-m1+1, n2-m2+1, n3-m3+1, nqf)
;qucold=0

;v1[*,*,130:*]=0.0
;v3[*,*,130:*]=0.0

; ei(where(ei lt 0.0   ))=3.0E+11
; ei(where(ei gt 1.0E14))=3.0E+11
; v1=(-7.0E+05) > v1 < 7.0E+05
; v2=(-7.0E+05) > v2 < 7.0E+05
; v3=(-7.0E+05) > v3 < 7.0E+05

; --- Spin-Up ---
;period_in_years=20.0
;rhd_manipulate_spinup, period_in_years*(365.2425*24.0*3600.0), m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3

; --- Generate new grid ---
;rhoavg=avg(avg(rho,0),0)
;dxb1=xb1[1]-xb1[0]
;dzb=replicate(0.5,n3-m3+1)
;;dzb[where(rhoavg gt 0.3e-4)]=1.00
;;dzb[where(rhoavg lt 0.5e-9)]=0.65
;dzb[where(alog10(rhoavg) gt -5.6)]=1.00
;dzb[where(alog10(rhoavg) lt -9.1)]=0.65
;;dzb[where(alog10(rhoavg) gt -3.5)]=1.00
;;dzb[where(alog10(rhoavg) lt -8.0)]=0.65
;dzb=dzb*dxb1
;help, ful, xc3new, P1, dzb
;rhd_manipulate_new1dgridgen, ful, xc3new, P1, dzb=dzb, j1max=560, aspect_min=0.5, aspect_max=1.0, dzchange=0.025, $
;   z1bot=xc3[0], z1top=xc3[-1], /input_type, z0in=xc3, dxbin=dxb1
;rhd_manipulate_interpolate, 0, 0, xc3new, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, indexstart='1', x3_command='newxc3'

;; --- Generate new grid ---
;;dxb1=(reform(ful.z.xb1))[1]-(reform(ful.z.xb1))[0]
;;dzb=replicate(0.50,n3-m3+1)
;;rhoavg=avg(avg(rho,0),0)
;;dzb[where(rhoavg gt 1e-4)]=1.0
;;dzb[where(rhoavg lt 1e-8)]=0.7
;;rhd_manipulate_new1dgridgen, ful, xc3new, P1, dzb=dzb, j1max=304, aspect_min=0.50, aspect_max=1.0, dzchange=0.03, $
;;   z1bot=(reform(ful.z.xc3))[0], z1top=(reform(ful.z.xc3))[n_elements(ful.z.xc3)-1]
;;rhd_manipulate_interpolate, 0, 0, xc3new, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, indexstart='1', x3_command='newxc3'


;dzb[184:*]=0.8
;dzb=dzb*dxb1
;rhd_manipulate_new1dgridgen, ful, xc3new, P1, dzb=dzb, j1max=280, aspect_min=0.40, aspect_max=1.0, dzchange=0.03, $
;   z1bot=(reform(ful.z.xc3))[0], z1top=(reform(ful.z.xc3))[n_elements(ful.z.xc3)-1]
;rhd_manipulate_interpolate, 0, 0, xc3new, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, indexstart='1', x3_command='newxc3'

;dxb1=(reform(ful.z.xb1))[1]-(reform(ful.z.xb1))[0]
;dzb=replicate(0.40,n_elements(ful.z.xc3))
;dzb[0:68]=1.0
;dzb[184:*]=0.8
;dzb=dzb*dxb1
;rhd_manipulate_new1dgridgen, ful, xc3new, P1, dzb=dzb, j1max=180, aspect_min=0.40, aspect_max=1.0, dzchange=0.03, $
;   z1bot=(reform(ful.z.xc3))[0], z1top=(reform(ful.z.xc3))[n_elements(ful.z.xc3)-1]
;rhd_manipulate_interpolate, 0, 0, xc3new, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, indexstart='1', x3_command='newxc3'
;delvar,ful

;xb3=xb3-xc3[0]
;xc3=xc3-xc3[0]

; dxb1=(reform(ful.z.xb1))(1)-(reform(ful.z.xb1))(0)
; dzb=replicate(0.5,n_elements(ful.z.xc3))
; dzb(0:50)=1.5
; dzb(89:*)=0.8
; dzb=dzb*dxb1
;rhd_manipulate_new1dgridgen, ful, xb3new, P1, dzb=dzb, j1max=121, aspect_min=0.28, aspect_max=1.7, dzchange=0.03, $
;   z1bot=(reform(ful.z.xb3))(0), z1top=(reform(ful.z.xb3))(n_elements(ful.z.xb3)-1)
;rhd_manipulate_interpolate, n1-m1+1, 0, xb3new, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='1'
;xc3=xc3-max(xb3)
;xb3=xb3-max(xb3)


;xfac=0.9
;xfac=0.8
;v1=v1*xfac
;v2=v2*xfac
;xc1=xc1*xfac & xb1=xb1*xfac
;xc2=xc2*xfac & xb2=xb2*xfac
;dxb1=(reform(xb1))(1)-(reform(xb1))(0)
;dzb=replicate(0.4,n_elements(xc3))
;dzb(0:60)=1.5
;dzb(0:20)=2.0
;dzb(111:*)=0.6
;dzb=dzb*dxb1
;rhd_manipulate_new1dgridgen, ful, xb3new, P1, dzb=dzb, j1max=171, aspect_min=0.17, aspect_max=2.0, dzchange=0.03, $
;  z1bot=(reform(xb3))(0), z1top=(reform(xb3))(n_elements(xb3)-1)
;rhd_manipulate_interpolate, 220, 220, xb3new, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='1'
;xc3=xc3-max(xb3)
;xb3=xb3-max(xb3)


;reffac=3
;dxb1=((reform(ful.z.xb1))(1)-(reform(ful.z.xb1))(0))/reffac
;n3new=54  ; 22, 90, 95, 30*reffac
;xb3new=dxb1*findgen(n3new) & xb3new=xb3new-max(xb3new)+xb3(n3)
;rhd_manipulate_interpolate, 0*(n1-m1+1)*reffac, 0*(n2-m2+1)*reffac, xb3new, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, indexstart='1'
;xc3=xc3-max(xb3)
;xb3=xb3-max(xb3)

;rhd_manipulate_interpolate, 200, 0, 0, m1,n1, m2,n2, m3,n3, xc1,xc2,xc3, xb1,xb2,xb3, rho, ei, v1,v2,v3, quc=quc, indexstart='box'


; --- Generate magnetic field arrays ---
;delvar, Bb1, Bb2, Bb3
;B3amp12=rebin(reform(sin(findgen(n1-m1+1)*(14.0*2.0*!pi/float(n1-m1+1))),n1-m1+1,1),n1-m1+1,n2-m2+1) * $
;        rebin(reform(sin(findgen(n2-m2+1)*(14.0*2.0*!pi/float(n2-m2+1))),1,n2-m2+1),n1-m1+1,n2-m2+1)
;  B3amp12=smooth2d(randomu(seed,n1-m1+1,n2-m2+1)*2.0, q_smooth=4,/periodic)
;  B3amp12=float(B3amp12-mean(double(B3amp12))) & B3amp12=B3amp12/sqrt(mean(B3amp12^2))
;;B3amp12=replicate(0.0, n1-m1+1, n2-m2+1) & B3amp12[(n1-m1+1)*1/4, (n2-m2+1)/2]=1.0 & B3amp12[(n1-m1+1)*3/4, (n2-m2+1)/2]=-1.0
;Bb1=replicate(0.0, n1-m1+2, n2-m2+1, n3-m3+1)
;Bb2=replicate(0.0, n1-m1+1, n2-m2+2, n3-m3+1)
;Bb3=replicate(0.0, n1-m1+1, n2-m2+1, n3-m3+2)
;Bb3=replicate(0.05/sqrt(4.0*!PI), n1-m1+1, n2-m2+1, n3-m3+2)
;Bb3=rebin(0.05/sqrt(4.0*!PI)*B3amp12, n1-m1+1, n2-m2+1, n3-m3+2)
;for j=40,50 do begin &$
;  Bb1[(n1-m1+1)/2-j  :(n1-m1+1)/2+j, (n2-m2+1)/2+j                , (n3-m3+1)/2-10:(n3-m3+1)/2+10]= 0.1 &$
;  Bb1[(n1-m1+1)/2-j  :(n1-m1+1)/2+j, (n2-m2+1)/2-j-1              , (n3-m3+1)/2-10:(n3-m3+1)/2+10]=-0.1 &$
;  Bb2[(n1-m1+1)/2+j                , (n2-m2+1)/2-j  :(n2-m2+1)/2+j, (n3-m3+1)/2-10:(n3-m3+1)/2+10]=-0.1 &$
;  Bb2[(n1-m1+1)/2-j-1              , (n2-m2+1)/2-j  :(n2-m2+1)/2+j, (n3-m3+1)/2-10:(n3-m3+1)/2+10]= 0.1 &$
;endfor ; i

; Initialise CHEM input
nquc=8
; Better to make 'name_species' and then add 'Number density of' as a prepend string after the call to 'co_initialabu()'
; This is just written out explicitly for now
name_quc='Number density of '+['H','H2','C','O','CO','CH','OH','metal']
name_species=['H','H2','C','O','CO','CH','OH','metal']  ; Same species as 'name_quc', passed to init abu stuff
unit_quc=replicate('cm^-3',nquc)
;quc=fltarr(szq[0],szq[1],szq[2],8)
quc=fltarr(n1-m1+1,n2-m2+1,n3-m3+1,8)

; metal abundance (as fraction of H density)
metal=0.1

; abundances at '$HOME/cobold/CHEM/abundance.dat'
; Make sure this corresponds with intended metallicity!
atomdata=co_loadatomdata()

; Calculate initial number density for each grid cell
for i3=0,n3-m3 do begin&$
  for i2=0,n2-m2 do begin&$
    for i1=0,n1-m1 do begin&$
    quc(i1,i2,i3,*)=float(co_initialabu(rho(i1,i2,i3),atomdata,name_species,metal=metal))&$
    endfor&$
  endfor&$
endfor&$

; quc=replicate(0.0, n1-m1+1, n2-m2+1, n3-m3+1, 1)
; name_quc=['Number density of CO']
; unit_quc=['1/cm^3']
;quc=replicate(0.0, n1-m1+1, n2-m2+1, n3-m3+1, 4)
;name_quc=['Dust moment 0', 'Dust moment 1', 'Dust moment 2', 'Dust moment 3']
;unit_quc=['1/cm^3', '1/cm^3', '1/cm^3', '1/cm^3']
;name_quc=['Dust moment 0', 'Dust moment 1']
;unit_quc=['1/cm^3', '1/cm^3']
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
