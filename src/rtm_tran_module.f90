module rtm_tran_module

!C-----------------------------------------------------------------------
!C !F90  rtm_tran_module                                                           
!C
!C !Description: 
!C   This module is to compute atmospheric transmittance.
!C
!C !Input  parameters
!C   none
!C 
!C !Output parameters
!C   none  
!C
!C  
!C !End
!C----------------------------------------------------------------------

! use modules
use constant
use names_module
use data_arrays_module, ONLY: zstd,pstd,tstd,wstd,ostd

!+++++++++++++++++++ step 1: define global variables +++++++++++++++++++
! |------|
! | none |
! |------|
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


contains
!+++++++++++++++++++ step 2: subroutines +++++++++++++++++++++++++++++++
subroutine tranvmodisd101(temp,wvmr,ozmr,theta,rco2,    &
                          craft,kban,jdet,              &
                          taut,iok)
! * MODIS band/detector 101-level fast transmittance routine
! .... version of 09.02.09

!        temp = profile of temperature ............ degK
!        wvmr = profile of H2O mixing ratio ....... g/kg
!        ozmr = profile of  O3 mixing ratio ....... ppmv
!       theta = local zenith angle ................ deg
!        rco2 = CO2 mixing ratio .................. ppmv
!       craft = TERRA, AQUA (upper or lower case)
!        kban = band number     (20...36)
!        jdet = detector number (0...10) [ Product Order ]
!                detector 0 is based on band-average response functions

!                taut = total transmittance (see note below)
!                 iok = 0 if successful, 1 if I/O problem

! * NOTE: for kban = 26, return-arrays are filled with 1.0

! * PLOD/PFAAST regression model based on LBLRTM line-by-line transmittances.
! * Input temperatures, and water-vapor and ozone mixing ratios, must
! *     be defined at the pressure levels in array 'pstd'
! *    (see block data 'reference_atmosphere').
! * Units: temperature, deg-K; water vapor, g/kg; ozone, ppmv.
! * Logical units 31-35 are used for coefficient files.
! * Component tau's are returned through common, product in 'taut'.

integer(kind=4) :: mets,kban
real(kind=4)    :: zena
integer(kind=4) :: l,k,i,j
integer(kind=4) :: lfac, nk, nl, nm, nr,nd,   &
                   nxc, ncc, lencc, lenccb,&
                   nxd, ncd, lencd, lencdb,&
                   nxo, nco, lenco, lencob,&
                   nxl, ncl, lencl, lenclb,&
                   nxs, ncs, lencs, lencsb,&
                   nxw
            
parameter (nd=10,nk=5,nl=101,nm=nl-1,koff=19,nr=17,lfac=1)
parameter (nxc= 4,ncc=nxc+1,lencc=ncc*nm,lenccb=lencc*lfac)
parameter (nxd= 8,ncd=nxd+1,lencd=ncd*nm,lencdb=lencd*lfac)
parameter (nxo= 9,nco=nxo+1,lenco=nco*nm,lencob=lenco*lfac)
parameter (nxl= 2,ncl=nxl+1,lencl=ncl*nm,lenclb=lencl*lfac)
parameter (nxs=11,ncs=nxs+1,lencs=ncs*nm,lencsb=lencs*lfac)
parameter (ndt=nd+1,nrps=nr*ndt,nxw=nxl+nxs)


      
!        common/stdatm/pstd(nl),tstd(nl),wstd(nl),ostd(nl)
!        common/taudwo/taud(nl),tauw(nl),tauo(nl)

real(kind=4), dimension(nl):: taud, tauw, tauo, taut
real(kind=4), dimension(nl) ::  temp, wvmr, ozmr

real(kind=4) :: coefd(ncd,nm,0:nd,nr),coefo(nco,nm,0:nd,nr),  &
                coefl(ncl,nm,0:nd,nr),coefs(ncs,nm,0:nd,nr),   &
                coefc(ncc,nm,0:nd,nr)
real(kind=4) :: bufc(lencc),bufd(lencd),bufo(lenco),  &
                bufl(lencl),bufs(lencs)
real(kind=4) :: pavg(nm),tref(nm),wref(nm),oref(nm)
real(kind=4) :: tavg(nm),wamt(nm),oamt(nm),secz(nm)
real(kind=4) :: tauc(nl),tlas(nl),wlas(nl),olas(nl)
real(kind=4) :: xdry(nxd,nm),xozo(nxo,nm),xwet(nxw,nm),xcon(nxc,nm)
        character*24 cfile(nk),xfile/'modisdet.com.101.xxx_end'/
        character*6 craft,cinit/'zzzzzz'/
        character*6 cbt/'TERRA'/,cba/'AQUA'/
        character*6 cst/'terra'/,csa/'aqua'/
        character*3 comp(nk)/'dry','ozo','wts','wtl','wco'/
        character*3 cbe/'big'/,cle/'lit'/
        character*20 path/'./data/plod/lit_end/'/
        integer*4 lengcf(nk)/lencdb,lencob,lencsb,lenclb,lenccb/
        integer*4 lengcx(nk)/lencd,lenco,lencs,lencl,lencc/
        integer*4 iuc(nk)
        logical big_endian,newang,newatm
        data tlas/nl*0./,wlas/nl*0./,olas/nl*0./,zlas/-999./
        
      !  secant(z)=1./cos(0.01745329*z)

     !   save
        if(craft.ne.cinit) then
           if(craft.eq.cbt.or.craft.eq.cst) then
              ksat=1
           elseif(craft.eq.cba.or.craft.eq.csa) then
              ksat=2
           else
              write(0,'(''In tran_modisd101 -- unknown spacecraft '',a6)') craft
              go to 200
           endif

! * determine which coefficient files to use
           !if(big_endian()) then
           !   xfile(18:20)=cbe
           !else
           !   xfile(18:20)=cle
           !endif
           xfile(18:20)=cle

! * define and open the files
           iux=30
           do m=1,nk
              iux=iux+1
              xfile(10:12)=comp(m)
              lencf=lengcf(m)
              
              !open(iux,file=path//xfile,recl=lencf,access='direct',status='old',err=200)
              open(iux,file=trim(code_root_path)//'coeff/plod/'//trim(xfile),recl=lencf,access='direct',status='old',err=200)
              iuc(m)=iux
              cfile(m)=xfile
           enddo

! * first read each file's fill-record for band 26/det 0
!               and verify satellite number stored in word 1
! * note: number of levels is in word 2, creation date (yyyyddd) is in word 3
           ikrec=nrps*(ksat-1)
           krecx=ikrec+7
           do k=1,nk
              lencx=lengcx(k)
              read(iuc(k),rec=krecx) (bufs(j),j=1,lencx)
              nsat=bufs(1)
              if(nsat.ne.ksat) then
                 xfile=cfile(k)
                 go to 100
              endif
           enddo

! * now read in the coefficients
           krec=ikrec
           do l=0,nd
              do k=1,nr
                  krec=krec+1
                  read(iuc(1),rec=krec) ((coefd(i,j,l,k),i=1,ncd),j=1,nm)
                  read(iuc(2),rec=krec) ((coefo(i,j,l,k),i=1,nco),j=1,nm)
                  read(iuc(3),rec=krec) ((coefs(i,j,l,k),i=1,ncs),j=1,nm)
                  read(iuc(4),rec=krec) ((coefl(i,j,l,k),i=1,ncl),j=1,nm)
                  read(iuc(5),rec=krec) ((coefc(i,j,l,k),i=1,ncc),j=1,nm)
              enddo
           enddo
           do k=1,nk
              close(iuc(k))
           enddo

           call conpir(pstd,tstd,wstd,ostd,nl,1,pavg,tref,wref,oref)
           cinit=craft
           iok=0
        endif

! * if ozone profile is null, put in std-atm
        if(ozmr(1).eq.0.) then
           do l=1,nl
              ozmr(l)=ostd(l)
           enddo
        endif

        dt=0.
        dw=0.
        do=0.
        do l=1,nl
           dt=dt+abs(temp(l)-tlas(l))
           tlas(l)=temp(l)
           dw=dw+abs(wvmr(l)-wlas(l))
           wlas(l)=wvmr(l)
           do=do+abs(ozmr(l)-olas(l))
           olas(l)=ozmr(l)
           taud(l)=1.0
           tauw(l)=1.0
           tauc(l)=1.0
           tauo(l)=1.0
           taut(l)=1.0
        enddo
        datm=dt+dw+do
! * if atmosphere has changed, convert atmospheric parameters
        newatm=datm.ne.0.
        if(newatm) then
           call conpir(pstd,temp,wvmr,ozmr,nl,1,pavg,tavg,wamt,oamt)
        endif

! * if angle has changed, recalculate secant
     !   newang=theta.ne.zlas
    !    if(newang) then
    !       zsec=secant(theta)
    !       do l=1,nm
    !          secz(l)=zsec
    !       enddo
    !       zlas=theta
    !    endif

if (theta /= zlas) then
   newang=.true.
endif

       if(newang) then
   !zsec=secant(zena) 1./cos(0.01745329*zena)
         zsec= 1./cos(0.01745329*theta)
         do l=1,nm
            secz(l)=zsec
         enddo
         zlas=zena
       endif

! * if atmosphere OR angle has changed, recalculate predictors
        if(newang.or.newatm) then
           call calpir(tref,wref,oref,tavg,wamt,oamt,pavg,secz,   &
                       nm,nxd,nxw,nxo,nxc,xdry,xwet,xozo,xcon)
        endif

        if(kban.eq.26) then
           do l=1,nl
              taud(l)=1.0
              tauo(l)=1.0
              tauw(l)=1.0
              taut(l)=1.0
           enddo
           return
        endif

        j=jdet
        k=kban-koff
! * dry
        call taudoc(ncd,nxd,nm,coefd(1,1,j,k),xdry,taud)
        if(rco2.ne.380.) then
! ... adjust for CO2 variation from model basis
           ratio=rco2/380.
           do l=1,nl
              if(taud(l).gt.0.0 .and. taud(l).lt.1.0) then
                 taud(l)=taud(l)**ratio
              endif
           enddo
        endif
! * ozo
        call taudoc(nco,nxo,nm,coefo(1,1,j,k),xozo,tauo)
! * wet
        call tauwtr(ncs,ncl,nxs,nxl,nxw,nm,coefs(1,1,j,k), &
                    coefl(1,1,j,k),xwet,tauw)
        call taudoc(ncc,nxc,nm,coefc(1,1,j,k),xcon,tauc)
        do l=1,nl
           tauw(l)=tauw(l)*tauc(l)
        enddo
! * total
        do l=1,nl
           taut(l)=taud(l)*tauo(l)*tauw(l)
        enddo
        return
100     write(0,'(''In tranvmodisd101 ... requested data for '', ''satellite '',i1/'' but read data for '', ''satellite '',i1,'' from file '',a24)') ksat,nsat,xfile
200     iok=1

        return
        
end subroutine tranvmodisd101
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine fy3dtrn101(temp,wvmr,ozmr,theta,kan,taut,iok)
! * Transmittance for FY-2E at 101 levels
! .... version of 18.11.05  (FY-2C)
! .... version of 15.09.11

! * LarrabeeStrow/HalWoolf/PaulVanDelst regression model based on
! *	LBLRTM line-by-line transmittances.
! * Input temperatures, and water-vapor and ozone mixing ratios, must
! *	be defined at the 101 pressure levels in array 'pstd' (see block
! *   data 'reference_atmosphere').
! * Units: temperature, deg-K; water vapor, g/kg; ozone, ppmv.
! * Logical unit numbers 71-75 are used for coefficient files.

! * Input
!	 temp = profile of temperature ........ degK
!	 wvmr = profile of H2O mixing ratio ... g/kg
!	 ozmr = profile of  O3 mixing ratio ... ppmv
!	theta = local zenith angle ............ deg
!c	  kan = channel number ................ 2 - 5
!	  kan = channel number ................ 1 - 4

! * Output
!	 taut = profile of total transmittance (components are returned through common)
!	    * = error return in case of coefficient-file I/O trouble

!	parameter (lfac=4,nk=5,nl=101,nm=nl-1,nr=5)
	parameter (lfac=4,nk=5,nl=101,nm=nl-1,nch=6,nr=nch+1)
	parameter (nxc= 4,ncc=nxc+1,lencc=ncc*nm,lenccb=lencc*lfac)
	parameter (nxd= 8,ncd=nxd+1,lencd=ncd*nm,lencdb=lencd*lfac)
	parameter (nxo= 9,nco=nxo+1,lenco=nco*nm,lencob=lenco*lfac)
	parameter (nxl= 2,ncl=nxl+1,lencl=ncl*nm,lenclb=lencl*lfac)
	parameter (nxs=11,ncs=nxs+1,lencs=ncs*nm,lencsb=lencs*lfac)
	parameter (nxw=nxl+nxs)
	common/stdatm/pstd(nl),tstd(nl),wstd(nl),ostd(nl)
	common/taudwo/taud(nl),tauw(nl),tauo(nl)
	dimension temp(*),wvmr(*),ozmr(*),taut(*)
	dimension coefd(ncd,nm,nr),coefo(nco,nm,nr),coefc(ncc,nm,nr)
	dimension coefl(ncl,nm,nr),coefs(ncs,nm,nr),iuc(nk)
	dimension pavg(nm),tref(nm),wref(nm),oref(nm)
	dimension tavg(nm),wamt(nm),oamt(nm),secz(nm)
	dimension tauc(nl),tlas(nl),wlas(nl),olas(nl)
	dimension xdry(nxd,nm),xozo(nxo,nm),xcon(nxc,nm),xwet(nxw,nm)
	character*7 cpath/'coeffs/'/
	character*14 cfile/'fy3dxxx101.dat'/
	character*3 comp(nk)/'dry','ozo','wco','wtl','wts'/
	integer*4 lencf(nk)/lencdb,lencob,lenccb,lenclb,lencsb/
	logical newang,newatm
	data init/1/,tlas/nl*0./,wlas/nl*0./,olas/nl*0./,zlas/-999./
	secant(z)=1./cos(0.01745329*z)

	if(init.ne.0) then
! * define and open the coefficient files
	   iux=92
	   do l=1,nk
	      cfile(5:7)=comp(l)
	      iux=iux+1
	      open(iux,file=trim(code_root_path)//'coeff/plod/'//cfile,recl=lencf(l)/4,access='direct',status='old',err=200)
	      iuc(l)=iux
	!     write(*,*) cpath//cfile
	   enddo
! * read in coefficients
	   do k=1,nr
	      krec=k
	      
	      read(iuc(1),rec=krec) ((coefd(i,j,k),i=1,ncd),j=1,nm)
	      read(iuc(2),rec=krec) ((coefo(i,j,k),i=1,nco),j=1,nm)
	      read(iuc(3),rec=krec) ((coefc(i,j,k),i=1,ncc),j=1,nm)
	      read(iuc(4),rec=krec) ((coefl(i,j,k),i=1,ncl),j=1,nm)
	      read(iuc(5),rec=krec) ((coefs(i,j,k),i=1,ncs),j=1,nm)
	   enddo
	   do l=1,nk
	      close(iuc(l))
	   enddo
	   
!	   write(16,*)coefd(:,:,1:9)
!	   write(17,*)coefo(:,:,1:9)
!	   write(18,*)coefc(:,:,1:9)
!	   write(19,*)coefl(:,:,1:9)
!	   write(20,*)coefs(:,:,1:9)
!	stop
!	    write(*,*) coefd(1,6,kan)
!	   write(*,*) coefd(4,60,kan),coefo(4,60,kan),coefc(4,60,kan)
!	   write(*,*) coefl(3,60,kan),coefs(4,60,kan)


! * initialize the reference profiles
	   call conpir(pstd,tstd,wstd,ostd,nl,1,pavg,tref,wref,oref)
	   init=0
	endif

	do j=1,nl
	   taud(j)=1.0
	   tauw(j)=1.0
	   tauc(j)=1.0
	   tauo(j)=1.0
	   taut(j)=1.0
	enddo
	

!	if(kan.eq.1) return

	dt=0.
	dw=0.
	do=0.
	do j=1,nl
	   dt=dt+abs(temp(j)-tlas(j))
	   tlas(j)=temp(j)
	   dw=dw+abs(wvmr(j)-wlas(j))
	   wlas(j)=wvmr(j)
	   do=do+abs(ozmr(j)-olas(j))
	   olas(j)=ozmr(j)
	enddo
	datm=dt+dw+do
	newatm=datm.ne.0.

	if(newatm) then
	   call conpir(pstd,temp,wvmr,ozmr,nl,1,pavg,tavg,wamt,oamt)
	endif

	newang=theta.ne.zlas
	if(newang) then
	   zsec=secant(theta)
	   do l=1,nm
	      secz(l)=zsec
	   enddo
	   zlas=theta
	endif

	if(newang.or.newatm) then
	   call calpir(tref,wref,oref,tavg,wamt,oamt,pavg,secz,    &
     		       nm,nxd,nxw,nxo,nxc,xdry,xwet,xozo,xcon)
	endif

	k=kan
	
! * dry
	call taudoc(ncd,nxd,nm,coefd(1,1,k),xdry,taud)
! * ozo
	call taudoc(nco,nxo,nm,coefo(1,1,k),xozo,tauo)
! * wet
! .. continuum
	call taudoc(ncc,nxc,nm,coefc(1,1,k),xcon,tauc)
	
! .. lines
	call tauwtr(ncs,ncl,nxs,nxl,nxw,nm,coefs(1,1,k),    &
                coefl(1,1,k),xwet,tauw)
	do j=1,nl
	   tauw(j)=tauw(j)*tauc(j)
	   
	enddo
! * total
	do j=1,nl
	   taut(j)=taud(j)*tauo(j)*tauw(j)
	  
	enddo
	return

!200	return
200 iok = 1
end subroutine fy3dtrn101
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module rtm_tran_module
