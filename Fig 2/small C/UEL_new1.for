       module kvisual
       implicit none
       real*8 sigout(10000,4,200)
C       public ::sigout
       save
       end module    

	   SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
		use kvisual
      INCLUDE 'ABA_PARAM.INC'
C
	  
       real(kind=8):: RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(NSVARS),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)



	   PARAMETER (maxDOF=2,ndof=2 , ndim=2 ,ninpt=4,nsdv=4,nsvint=200,ntens=4)!Here nsdv only shows the number of stress
       real(kind=8):: dshape(ndim,4),xjac(maxDOF,maxDOF),xjaci(maxDOF,maxDOF)
       real(kind=8):: kmat(2*NNODE,2*NNODE),C(3,3),STRAN(3),shapef(4)
	   real(kind=8):: bmata(3,2*NNODE),bmataT(2*NNODE,3)
       real(kind=8):: deps(3),stress(3),statevLocal(200)!All the values are stored
	   real(kind=8):: bmatb(8,16),bmatbT(8,16),bmatc(4,8),bmatcT(8,4)
	   real(kind=8):: nmatb(16,4),nmatbT(4,16)
	   real(kind=8):: XYDerivation(ndim,4),force(40),DDSDDE(3,3),ff(40)
	   real(kind=8) :: fd(8,1)
       real(kind=8) :: feta(16,1)
       real(kind=8) :: flambda(16,1)
       real(kind=8) :: K_f(40, 40)
	   real(kind=8) :: real_u(8),real_du(8,1)
       real(KIND=8) :: Mazarslocal(190),Sprainlocal(2)
	   real(KIND=8) :: sp_energy,st_energy
      
C       WRITE(6,*), "JELEM=", JELEM
c       WRITE(6,*), "KINC=", KINC
C       WRITE(6,*), "U32=", U
C       WRITE(6,*), "DU33=", DU(1:40,1)
C       integer :: NTENS = 3
		rhs(:, 1)=0.d0
		AMATRX=0.d0
		C=0.d0		
C		write(6,*), "PROPS=",PROPS(1:6)
c        write(6,*), "NSVARS=",NSVARS	
c	  write(*,*) 'DDSDDE is ', kind(DDSDDE)	
        force=0.d0
	    DDSDDE = 0.d0
        Mazarslocal = 0.0d0
        !Elasticity matrix
c		E = props(1)
c		xnu = props(2)
c		xa = E/(1.d0- xnu**2)
c		xb = xnu
c		xc = (1.d0-xnu)/2.d0
		
c        DTIME = 0.001
        do kintk = 1, ninpt	 
C		write(6,*), "kintk=",kintk	
C		write(6,*), "time=",time(1)	
		deps=0.	
		statevLocal=0.0
        call shapeFunc(kintk,ninpt,NNODE,ndim,shapef,dshape)
		call shapeFunc2(shapef,nmatb,nmatbT)
C		WRITE(6,*), "shapef=", shapef
C		WRITE(6,*), "dshape=", dshape
c		WRITE(6,*), "COORDS=", COORDS
C		pause
	    djac = 1.d0
        call jacobian(ndim,COORDS,NNODE,djac,dshape,xjaci)
        	    
        call Bmatrix(xjaci,dshape,NNODE,ndim,bmata,bmataT)

        call Bmatrix2(xjaci,dshape,NNODE,ndim,bmatb,bmatbT,bmatc,bmatcT)

c        write(*,*) 'force is ', kind(force)	
c         write(6,*), "shapef=",shapef
C         write(6,*), "nmatb=",nmatb
C         write(6,*), "bmata=",bmata
C		 write(6,*), "bmatb=",bmatb
C		 write(6,*), "bmatc=",bmatc
c		 PAUSE
c        write(*,*),"bmat_row",shape(bmat, 1)
c        write(*,*),"bmat_column",shape(bmat, 2)
	    call statevar(kintk,nsvint,SVARS,NSVARS,statevLocal,1)
c          write(*,*) 'stress is ', kind(stress)	
		stress=0.d0
        STRAN=0.d0
		!pass the statevLocal values to the current calculation
		stress(1:3) = statevLocal(1:3)
		STRAN(1:3) = statevLocal(4:6)
        Mazarslocal(1:190) = statevLocal(7:196)  
		Sprainlocal(1:2) = statevLocal(197:198)  
		sp_energy = statevLocal(199)
		st_energy = statevLocal(200)  
c        write(*,*) 'statevLocal(8:197)  is ', kind(statevLocal)	
c	     write(*,*) 'Mazarslocal is ', kind(Mazarslocal)		
c        write(6,*), "stress80=",stress
c        write(6,*), "Mazarslocal=",Mazarslocal

        real_du = 0d0
        real_du(1,1) = DU(1,1)
		real_du(2,1) = DU(2,1)
		real_du(3,1) = DU(11,1)
		real_du(4,1) = DU(12,1)
		real_du(5,1) = DU(21,1)
		real_du(6,1) = DU(22,1)
		real_du(7,1) = DU(31,1)
		real_du(8,1) = DU(32,1)

		call straininc(ntens,ndof,ndim,NNODE,8,bmata,real_du,deps)	!MLVARX=8
c		write(6,*), "U105=",U(:)	
c        write(6,*), "DU83=",DU(:,1)		
c        STRAN = STRAN + deps
C        write(6,*), "real_du=",real_du
C        write(6,*), "STRAN old=",STRAN	
C        write(6,*), "deps85=",deps	
c        write(*,*), "bmat87=",bmat	
c        write(6,*), "UELMazarslocal before=",Mazarslocal	

c        IF (JELEM == 101 ) THEN
c        PROPS(4) = PROPS(4)*0.80d0
c	    END IF
c		STRAN = STRAN + deps
c           call M7FMATERIAL(stress,M7local,DDSDDE,STRAN, deps,
c     $       TIME,DTIME, 3,190,PROPS,NPROPS,COORDS)
C	    IF (JELEM == 6 ) THEN
C        call M7FMATERIAL2(stress,Mazarslocal,DDSDDE,STRAN, deps,
C     $       TIME,DTIME, 3,2,PROPS,NPROPS,COORDS) !use M7 laws NTENS=3, STATE=190
C	    else
        call M7FMATERIAL(stress,Mazarslocal,DDSDDE,STRAN, deps,
     $       TIME,DTIME, 3,2,PROPS,NPROPS,COORDS) !use M7 laws NTENS=3, STATE=190
C		END IF
C       write(6,*), "UELMazarslocal after=",Mazarslocal	
c     !NEED STRAIN TO CALCULATE
c        IF (JELEM == 101) THEN
c	 	call M7FMATERIAL2(stress,Mazarslocal,DDSDDE,STRAN, deps,
c     $       TIME,DTIME, 3,190,PROPS,NPROPS,COORDS) !use M7 laws NTENS=3, STATE=190
C        stress = stress * 0.9D0
c	    END IF
        STRAN = STRAN + deps
c        call Bmatrix(xjaci,dshape,NNODE,ndim,bmat,bmatT)

C		stress=stress+matmul(DDSDDE,deps) !non linear shoud not use this LINE

c        write(6,*), "stress99=",stress
		 call SprainE(Sprainlocal,PROPS,NPROPS,bmata,
	1        bmatb,bmatc,nmatb,NDOFEL,
	2        MLVARX,U,DU(:,1),DDSDDE,shapef,
	3       stress,strain,fd,feta,flambda,K_f,COORDS,sp_energy)
		IF (JELEM == 491) THEN
C		write(6,*), "sp_energy=",sp_energy 
C	 	call Sprain2(bmata,bmatb,bmatc,nmatb,NDOFEL,MLVARX,U,DU(:,1),DDSDDE,shapef,stress,fd,feta,flambda,K_f)
	    END IF
		force=0.d0
C        stress = RESHAPE(stress, [3, 1])
        force(1:8)=fd(1:8,1)
C		write(6,*), "fd=", fd(:,1)
		force(9:24)=feta(1:16,1)
C		write(6,*), "feta=", feta(:,1)
		force(25:40)=flambda(1:16,1)
C		write(6,*), "flambda=", flambda(:,1)
c		force=matmul(bmataT,stress)
		force=force*djac
        call finalforce(force,ff)
c        write(6,*), "force100=",force	
		do k1=1,40
				RHS(k1, 1) = RHS(k1, 1) - ff(k1)
		end do
c		write(6,*), "RHS104=", RHS(:,1)
        


        AMATRX = AMATRX + K_f*djac
C		write(6,*), "djac=",djac
C		write(6,*), "AMATRX=",AMATRX
c		pause
c        write(6,*), "SIG=",SIG	
		st_energy =0.0d0
        DO i =1,3
		st_energy=st_energy+0.5D0*STRAN(i)*stress(i)
		END DO   
        statevLocal(1:3) = stress(1:3)
        statevLocal(4:6) = STRAN(1:3)
c        write(*,*) 'statevLocal is ', kind(statevLocal)
        statevLocal(7:196) = Mazarslocal(1:190)  
		statevLocal(197:198) = Sprainlocal(1:2)         
c        write(6,*), "statevLocal=",statevLocal
        statevLocal(199) = sp_energy 
c		if (st_energy /= 0.0d0) then
c		statevLocal(200) = sp_energy/st_energy 
c		else
c		statevLocal(200) = 0.0d0
c		end if
		statevLocal(200) = st_energy 
        call statevar(kintk,nsvint,SVARS,NSVARS,statevLocal,0)
c        write(6,*), "kintk=",kintk
c        write(6,*), "SVARS=",SVARS
	

		do k1=1,3
			sigout(jelem,kintk,k1)=stress(k1)
		end do
        do k1=1,3
			sigout(jelem,kintk,3+k1)=STRAN(k1)
		end do
		do k1=1,190
			sigout(jelem,kintk,6+k1)=Mazarslocal(k1)
		end do
		do k1=1,2
			sigout(jelem,kintk,196+k1)=Sprainlocal(k1)
		end do		
		sigout(jelem,kintk,199) = sp_energy 
		sigout(jelem,kintk,200) = sp_energy/(st_energy+ sp_energy)
		
		end do       

c        pause
      RETURN
      END
	  
********************************************************************
      subroutine finalforce(force, ff)
      real(kind=8):: force(40),ff(40)
      integer(kind=8):: ind(40)
      ind =   (/1,2,9,10,11,12,25,26,27,28,
     $          3,4,13,14,15,16,29,30,31,32,
     $          5,6,17,18,19,20,33,34,35,36,
     $          7,8,21,22,23,24,37,38,39,40/)
      do i = 1,40
          ff(i)=force(ind(i))
      end do
      
      end subroutine finalforce	  
	  
******************************************************
** shape functions                                   *
******************************************************	
		subroutine shapeFunc(kintk,ninpt,NNODE,ndim,shapef,dshape)
		
		include 'aba_param.inc'
		
		parameter (gaussCoord=0.577350269189626)
		real(kind=8):: shapef(4),dshape(ndim,4),coord24(2,4),r,s
		
		data coord24  /-1. , -1. ,
	2				  -1. , 1. ,
	3				   1. , -1.,
	4				   1. , 1./
			
		shapef=0.
		dshape=0.
C		write(6,*), "coord24=",coord24
		r=coord24(1,kintk)*gaussCoord
		s=coord24(2,kintk)*gaussCoord
		
		! shape functions
		
		shapef(1)=0.25d0*(1.-r)*(1.-s)
		shapef(2)=0.25d0*(1.+r)*(1.-s)
        shapef(3)=0.25d0*(1.+r)*(1.+s)
        shapef(4)=0.25d0*(1-r)*(1.+s)
	
		! derivation to r
		
		dshape(1,1)=0.25d0*(-1.+s)
		dshape(1,2)=0.25d0*(1.-s)
        dshape(1,3)=0.25d0*(1.+s)
        dshape(1,4)=0.25d0*(-1.-s)	
	
	
	    ! derivation to s
		
		dshape(2,1)=0.25d0*(-1.+r)
		dshape(2,2)=0.25d0*(-1.-r)
        dshape(2,3)=0.25d0*(1.+r)
        dshape(2,4)=0.25d0*(1.-r)
	
	
		return
		end
******************************************************
** shape functions                                   *
******************************************************	
		subroutine shapeFunc2(shapef,nmatb,nmatbT)
		
c		include 'aba_param.inc'
		
		real(kind=8):: shapef(4),nmatb(4,16),nmatbT(16,4)
		nmatb = 0.d0
		do I = 1,4
		nmatb(1,1+4*(I-1))=shapef(I)
		nmatb(2,2+4*(I-1))=shapef(I)
		nmatb(3,3+4*(I-1))=shapef(I)
		nmatb(4,4+4*(I-1))=shapef(I)
		end do	
C		write(6,*), "nmatb in shapeFunc2=",nmatb
		do i=1,4
			do j=1,16
				nmatbT(j,i)=nmatb(i,j)
			end do
		end do
	
		return
		end
		
******************************************************
** Jacobian                                          *
******************************************************
		subroutine jacobian(ndim,coords,NNODE,djac,dshape,xjaci)
		
		include 'aba_param.inc'
		
		real(kind=8):: xjac(ndim,ndim),xjaci(ndim,ndim)
		real(kind=8):: coords(2,NNODE),dshape(ndim,NNODE)

		xjac=0.d0
		xjaci=0.d0

		
        do inod= 1, NNODE
          do kdim = 1, ndim
            do jdim = 1, ndim
              xjac(jdim,kdim) = xjac(jdim,kdim) + 
     1        dshape(jdim,inod)*coords(kdim,inod)      
            end do
          end do 
        end do
		
		
		djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1)
		 
		if (djac .gt. 0.) then
			xjaci(1,1)=xjac(2,2)/djac
			xjaci(2,2)=xjac(1,1)/djac
			xjaci(1,2)=-xjac(1,2)/djac
			xjaci(2,1)=-xjac(2,1)/djac			
		end if
		
		return
		end
******************************************************
** B-matrix                                          *
******************************************************
		subroutine Bmatrix(xjaci,dshape,NNODE,ndim,bmat,bmatT)
		
		include 'aba_param.inc'
		
		real(kind=8):: xjaci(ndim,ndim),dshape(ndim,4),XYDerivation(ndim,4)
		real(kind=8):: bmat(3,2*NNODE),bmatT(2*NNODE,3)

		XYDerivation=0.d0

		bmat=0.d0
		bmatT=0.d0


		XYDerivation=matmul(xjaci,dshape)

				

		bmat(1,1)=XYDerivation(1,1)
		bmat(1,3)=XYDerivation(1,2)
		bmat(1,5)=XYDerivation(1,3)
		bmat(1,7)=XYDerivation(1,4)
		
		bmat(2,2)=XYDerivation(2,1)
		bmat(2,4)=XYDerivation(2,2)
		bmat(2,6)=XYDerivation(2,3)
		bmat(2,8)=XYDerivation(2,4)
		
		bmat(3,1)=XYDerivation(2,1)
		bmat(3,2)=XYDerivation(1,1)
		bmat(3,3)=XYDerivation(2,2)
		bmat(3,4)=XYDerivation(1,2)
		bmat(3,5)=XYDerivation(2,3)
		bmat(3,6)=XYDerivation(1,3)
		bmat(3,7)=XYDerivation(2,4)
		bmat(3,8)=XYDerivation(1,4)
		
		do i=1,3
			do j=1,2*NNODE
				bmatT(j,i)=bmat(i,j)
			end do
		end do
		
		return
		end
******************************************************
** B-matrix                                          *
******************************************************
		subroutine Bmatrix2(xjaci,dshape,NNODE,ndim,bmatb,bmatbT,bmatc,bmatcT)
		
		include 'aba_param.inc'
		
		real(kind=8):: xjaci(ndim,ndim),dshape(ndim,4),XYDerivation(ndim,4)
		real(kind=8):: bmatb(8,16),bmatbT(16,8),bmatc(4,8),bmatcT(8,4)

		XYDerivation=0.d0

		bmatb=0.d0
		bmatbT=0.d0
        bmatc=0.d0
		bmatcT=0.d0

		XYDerivation=matmul(xjaci,dshape)
C        WRITE(6,*),"XYDerivation",XYDerivation 
				
        do I = 1,4
		bmatb(1,1+4*(I-1))=XYDerivation(1,I)
		bmatb(2,1+4*(I-1))=XYDerivation(2,I)
		bmatb(3,2+4*(I-1))=XYDerivation(1,I)
		bmatb(4,2+4*(I-1))=XYDerivation(2,I)
		bmatb(5,3+4*(I-1))=XYDerivation(1,I)
		bmatb(6,3+4*(I-1))=XYDerivation(2,I)
		bmatb(7,4+4*(I-1))=XYDerivation(1,I)
		bmatb(8,4+4*(I-1))=XYDerivation(2,I)
		end do	
			
		do i=1,8
			do j=1,16
				bmatbT(j,i)=bmatb(i,j)
			end do
		end do

        do I = 1,4
		bmatc(1,1+2*(I-1))=XYDerivation(1,I)
		bmatc(2,1+2*(I-1))=XYDerivation(2,I)
		bmatc(3,2+2*(I-1))=XYDerivation(1,I)
		bmatc(4,2+2*(I-1))=XYDerivation(2,I)
		end do	
			
		do i=1,4
			do j=1,8
				bmatcT(j,i)=bmatc(i,j)
			end do
		end do

		
		return
		end
		
c*****************************************************************
		  subroutine straininc(ntens,ndof,ndim,NNODE,mlvarx,bbmat,du,deps)

		  !include 'aba_param.inc'
		  INTEGER:: ndof,NNODE,mlvarx
		  real(kind=8):: deps(3),bbmat(3,2*NNODE),du(mlvarx,1),xdu(2) 
		  
			
		  deps = 0.d0
		  xdu=0.d0
			
		  do nodi = 1, NNODE
			   
		   incr_row = (nodi - 1)*ndof
			
		   do i = 1, ndof         
				xdu(i)= du(i + incr_row,1)
		   end do
			dNidx = bbmat(1,1 + incr_row)
			dNidy = bbmat(2,1 + incr_row + 1)
			
		   deps(1) = deps(1) + dNidx*xdu(1)
		   deps(2) = deps(2) + dNidy*xdu(2)
		   deps(3) = deps(3) + dNidy*xdu(1) + dNidx*xdu(2)

          
C		  WRITE(*,*),"dNidx",dNidx
C          WRITE(*,*),"dNidy",dNidy	

		  end do
c          WRITE(*,*),"bmat299",bbmat 
		  return
		  end


		
************************************************************************		
		 subroutine statevar(npt,nsvint,statev,NSVARS,statev_ip,icopy)
	

		  include 'aba_param.inc'

		  real(kind=8):: statev(NSVARS),statev_ip(nsvint)

		  isvinc= (npt-1)*nsvint     

		  if (icopy .eq. 1) then

			do i = 1, nsvint
			  statev_ip(i)=statev(i+isvinc)
			end do
c
		  else
c
			do i = 1, nsvint
			  statev(i+isvinc)=statev_ip(i)
			end do
		  end if

		  return
		  end
************************************************************************
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
c 
c     NOTE: developed by Hoang T Nguyen Aug 2020
c
c 
C-----  Use single precision on Cray by
C     (1) deLeting the statement "IMPLICIT*8 (A-H,O-Z)";
C     (2) changing "REAL*8 FUNCTION" to "FUNCTION";
C     (3) changing DOuble precision functions DSIGN to SIGN.
C
C-----  Subroutines:
C        use kvisual
      INCLUDE 'aba_param.inc'
c
      CHARACTER*80 CMNAME

      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      IF(CMNAME(1:11).EQ."DECOY") THEN
         CALL DECOY(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      END IF

      RETURN
      END SUBROUTINE UMAT
C **********************************************************************
C *** SUBROUTINE DECOY     *********************************************
C **********************************************************************
       subroutine DECOY(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
C
	  use kvisual
      include 'aba_param.inc'
C
      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),
     2 ddsddt(ntens),drplde(ntens),
     3 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     4 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
	  PARAMETER (maxDOF=2 , ndof=2 , ndim=2 ,ninpt=4,nelem=1940)
	  integer kelem
C      REAL*8, INTENT(IN) :: sigout_arg(10000,4,4)


		ddsdde=0.d0

      kelem=int(noel-nelem)
	 
      do k1 = 1, 6
        statev(k1) = sigout(kelem,npt,k1)
c        write(6,*), "k1=",k1
      end do 
      do k1 = 7, 8
        statev(k1) = sigout(kelem,npt,190+k1)
c        write(6,*), "k1=",k1
      end do 
	  statev(9)  =  sigout(kelem,npt,199)
      statev(10) =  sigout(kelem,npt,200)
      return
      end subroutine DECOY
c      include 'Mazars_Damage_Model.for'
      include 'M7fMATERIAL.for'
C	  include 'M7fMATERIAL2.for'
C	  include 'Sprain.for'
	  include 'SprainE.for'
c	  include 'Bilinear_Damage_Model.for'
c	  include 'Bilinear_Damage_Model_weak.for'