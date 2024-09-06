
C ADJUSTING PARAMETERS TO FIT GIVEN MATERIAL TEST DATA IN THE 
C INPUTPARAMS() SUBROUTINE
C
C This is the implicit user-defined subroutine for M7 version of microplane
C model for concrete

C **********************************************************************
C *** SUBROUTINE SETSYSTEM *********************************************
C **********************************************************************
 
      SUBROUTINE SETSYSTEM()
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL(KIND=8), PARAMETER :: PI=3.1415926535897932384626433832795d0
      INTEGER, PARAMETER ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      COMMON /KCOMMON_Vars/young,poisson,k_1,th_del,unit_conv
      REAL(KIND=8) :: young,poisson,k_1,th_del,unit_conv,zeta,zeta0
      COMMON /KM7f_Microplane_1/ qn, ql, qm, w 
      REAL(KIND=8), DIMENSION(1:6,1:np) :: qn, ql, qm
      REAL(KIND=8), DIMENSION(1:np)     :: w
      SAVE    :: /KM7f_Microplane_1/
      SAVE    :: /KCOMMON_Vars/
      INTEGER:: jp, ij(1:2,1:6), i, j, k, allocstat
      INTEGER:: rs_size
      REAL(KIND=8), DIMENSION(1:4,1:np) :: te
      REAL(KIND=8), DIMENSION(1:3) :: xn, xm, xl, rand_vec 
      REAL(KIND=8) :: lengthn, lengthm, lengthl
C
          qn = 0.d0
          ql = 0.d0
          qm = 0.d0
          w = 0.d0
          te = 0.d0
          
C     GENERATE THE TABLE TO NUMERICALLY CALCULATE THE SPHERICAL INTEGRAL
C     This depends on the number of integration points
      ij=RESHAPE((/1,1,2,2,3,3,1,2,1,3,2,3/),(/2,6/))
      IF (np.eq.21) THEN
      te = RESHAPE(
     $(/0.000000000000D+00,0.000000000000D+00,1.000000000000D+00,
     $        2.652142440930D-02,
     $ 0.000000000000D+00,1.000000000000D+00,0.000000000000D+00,
     $        2.652142440930D-02,
     $ 1.000000000000D+00,0.000000000000D+00,0.000000000000D+00,
     $        2.652142440930D-02,
     $ 0.000000000000D+00,7.071067811870D-01,7.071067811870D-01,
     $        1.993014763120D-02,
     $ 0.000000000000D+00,-7.071067811870D-01,7.071067811870D-01,
     $        1.993014763120D-02,
     $ 7.071067811870D-01,0.000000000000D+00,7.071067811870D-01,
     $        1.993014763120D-02,
     $-7.071067811870D-01,0.000000000000D+00,7.071067811870D-01,
     $        1.993014763120D-02,
     $ 7.071067811870D-01,7.071067811870D-01,0.000000000000D+00,
     $        1.993014763120D-02,
     $-7.071067811870D-01,7.071067811870D-01,0.000000000000D+00,
     $        1.993014763120D-02,
     $ 8.360955967490D-01,3.879073040670D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $-8.360955967490D-01,3.879073040670D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $ 8.360955967490D-01,-3.879073040670D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $-8.360955967490D-01,-3.879073040670D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $ 3.879073040670D-01,8.360955967490D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $-3.879073040670D-01,8.360955967490D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $ 3.879073040670D-01,-8.360955967490D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $-3.879073040670D-01,-8.360955967490D-01,3.879073040670D-01,
     $        2.507123674870D-02,
     $ 3.879073040670D-01,3.879073040670D-01,8.360955967490D-01,
     $        2.507123674870D-02,
     $-3.879073040670D-01,3.879073040670D-01,8.360955967490D-01,
     $        2.507123674870D-02,
     $ 3.879073040670D-01,-3.879073040670D-01,8.360955967490D-01,
     $        2.507123674870D-02,
     $-3.879073040670D-01,-3.879073040670D-01,8.360955967490D-01,
     $        2.507123674870D-02/)
     $,(/4,np/))
      ELSE IF (np.eq.28) THEN

      te = RESHAPE(
     $(/.577350259D0, .577350259D0, .577350259D0, .0160714276D0,
     $.577350259D0, .577350259D0,-.577350259D0, .0160714276D0,
     $.577350259D0,-.577350259D0, .577350259D0, .0160714276D0,
     $.577350259D0,-.577350259D0,-.577350259D0, .0160714276D0,
     $.935113132D0, .250562787D0, .250562787D0, .0204744730D0,
     $.935113132D0, .250562787D0,-.250562787D0, .0204744730D0,
     $.935113132D0,-.250562787D0, .250562787D0, .0204744730D0,
     $.935113132D0,-.250562787D0,-.250562787D0, .0204744730D0,
     $.250562787D0, .935113132D0, .250562787D0, .0204744730D0,
     $.250562787D0, .935113132D0,-.250562787D0, .0204744730D0,
     $.250562787D0,-.935113132D0, .250562787D0, .0204744730D0,
     $.250562787D0,-.935113132D0,-.250562787D0, .0204744730D0,
     $.250562787D0, .250562787D0, .935113132D0, .0204744730D0,
     $.250562787D0, .250562787D0,-.935113132D0, .0204744730D0,
     $.250562787D0,-.250562787D0, .935113132D0, .0204744730D0,
     $.250562787D0,-.250562787D0,-.935113132D0, .0204744730D0,
     $.186156720D0, .694746614D0, .694746614D0, .0158350505D0,
     $.186156720D0, .694746614D0,-.694746614D0, .0158350505D0,
     $.186156720D0,-.694746614D0, .694746614D0, .0158350505D0,
     $.186156720D0,-.694746614D0,-.694746614D0, .0158350505D0,
     $.694746614D0, .186156720D0, .694746614D0, .0158350505D0,
     $.694746614D0, .186156720D0,-.694746614D0, .0158350505D0,
     $.694746614D0,-.186156720D0, .694746614D0, .0158350505D0,
     $.694746614D0,-.186156720D0,-.694746614D0, .0158350505D0,
     $.694746614D0, .694746614D0, .186156720D0, .0158350505D0,
     $.694746614D0, .694746614D0,-.186156720D0, .0158350505D0,
     $.694746614D0,-.694746614D0, .186156720D0, .0158350505D0,
     $.694746614D0,-.694746614D0,-.186156720D0, .0158350505D0/)
     $,(/4,np/))     
      ELSE IF (np.eq.37) THEN
      te = RESHAPE(
     $(/0.000000000000D+00,0.000000000000D+00,1.000000000000D+00,
     $        1.072388573030D-02,
     $ 0.000000000000D+00,1.000000000000D+00,0.000000000000D+00,
     $        1.072388573030D-02,
     $ 1.000000000000D+00,0.000000000000D+00,0.000000000000D+00,
     $        1.072388573030D-02,
     $ 0.000000000000D+00,7.071067811870D-01,7.071067811870D-01,
     $        2.114160951980D-02,
     $ 0.000000000000D+00,-7.071067811870D-01,7.071067811870D-01,
     $        2.114160951980D-02,
     $ 7.071067811870D-01,0.000000000000D+00,7.071067811870D-01,
     $        2.114160951980D-02,
     $-7.071067811870D-01,0.000000000000D+00,7.071067811870D-01,
     $        2.114160951980D-02,
     $ 7.071067811870D-01,7.071067811870D-01,0.000000000000D+00,
     $        2.114160951980D-02,
     $-7.071067811870D-01,7.071067811870D-01,0.000000000000D+00,
     $        2.114160951980D-02,
     $ 0.000000000000D+00,3.089512677750D-01,9.510778696510D-01,
     $        5.355055908370D-03,
     $ 0.000000000000D+00,-3.089512677750D-01,9.510778696510D-01,
     $        5.355055908370D-03,
     $ 0.000000000000D+00,9.510778696510D-01,3.089512677750D-01,
     $        5.355055908370D-03,
     $ 0.000000000000D+00,-9.510778696510D-01,3.089512677750D-01,
     $        5.355055908370D-03,
     $ 3.089512677750D-01,0.000000000000D+00,9.510778696510D-01,
     $        5.355055908370D-03,
     $-3.089512677750D-01,0.000000000000D+00,9.510778696510D-01,
     $        5.355055908370D-03,
     $ 9.510778696510D-01,0.000000000000D+00,3.089512677750D-01,
     $        5.355055908370D-03,
     $-9.510778696510D-01,0.000000000000D+00,3.089512677750D-01,
     $        5.355055908370D-03,
     $ 3.089512677750D-01,9.510778696510D-01,0.000000000000D+00,
     $        5.355055908370D-03,
     $-3.089512677750D-01,9.510778696510D-01,0.000000000000D+00,
     $        5.355055908370D-03,
     $ 9.510778696510D-01,3.089512677750D-01,0.000000000000D+00,
     $        5.355055908370D-03,
     $-9.510778696510D-01,3.089512677750D-01,0.000000000000D+00,
     $        5.355055908370D-03,
     $ 8.805355183100D-01,3.351545919390D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $-8.805355183100D-01,3.351545919390D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $ 8.805355183100D-01,-3.351545919390D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $-8.805355183100D-01,-3.351545919390D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $ 3.351545919390D-01,8.805355183100D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $-3.351545919390D-01,8.805355183100D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $ 3.351545919390D-01,-8.805355183100D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $-3.351545919390D-01,-8.805355183100D-01,3.351545919390D-01,
     $        1.677709091560D-02,
     $ 3.351545919390D-01,3.351545919390D-01,8.805355183100D-01,
     $        1.677709091560D-02,
     $-3.351545919390D-01,3.351545919390D-01,8.805355183100D-01,
     $        1.677709091560D-02,
     $ 3.351545919390D-01,-3.351545919390D-01,8.805355183100D-01,
     $        1.677709091560D-02,
     $-3.351545919390D-01,-3.351545919390D-01,8.805355183100D-01,
     $        1.677709091560D-02,
     $ 5.773502691900D-01,5.773502691900D-01,5.773502691900D-01,
     $        1.884823095080D-02,
     $-5.773502691900D-01,5.773502691900D-01,5.773502691900D-01,
     $        1.884823095080D-02,
     $ 5.773502691900D-01,-5.773502691900D-01,5.773502691900D-01,
     $        1.884823095080D-02,
     $-5.773502691900D-01,-5.773502691900D-01,5.773502691900D-01,
     $        1.884823095080D-02/)
     $,(/4,np/))

      END IF   
C--------------------------------------------
C     Assemble tensors from direction cosines
C--------------------------------------------
      qn = 0.d0
      ql = 0.d0
      qm = 0.d0
      w = 0.d0
      CALL random_seed(size=rs_size)
      DO jp=1,np 
          w(jp) = te(4,jp)*6.0D0 
          xn(1) = te(3,jp) 
          xn(2) = te(2,jp) 
          xn(3) = te(1,jp) 
C     Pick up a ranDOm vector that is perpendicular to xn
          lengthm = 0.d0
          DO WHILE (lengthm .lt. epsilon(lengthm))
              CALL random_number(rand_vec)
              xm = rand_vec - dot_product(xn,rand_vec)*xn
              lengthm = sqrt(dot_product(xm,xm))
          END DO
          xm = xm/lengthm
C     From xn and xm calculate xl from cross product
          xl(1) = xn(2)*xm(3)-xn(3)*xm(2) 
          xl(2) = xn(3)*xm(1)-xn(1)*xm(3) 
          xl(3) = xn(1)*xm(2)-xn(2)*xm(1) 
          lengthl = sqrt(dot_product(xl,xl))
          xl=xl/lengthl
          lengthn = sqrt(dot_product(xn,xn))
          lengthm = sqrt(dot_product(xm,xm))
          lengthl = sqrt(dot_product(xl,xl))
          DO k=1,6 
              i=ij(1,k) 
              j=ij(2,k) 
              qn(k,jp) = xn(i)*xn(j) 
              qm(k,jp) = 0.5D0*(xn(i)*xm(j)+xn(j)*xm(i))
              ql(k,jp) = 0.5D0*(xn(i)*xl(j)+xn(j)*xl(i)) 
          END DO 
      END DO 
      RETURN
      END SUBROUTINE SETSYSTEM 

C *******************************************************************
C *** SUBROUTINE M7fMATERIAL *****************************************
C *******************************************************************

      SUBROUTINE M7FMATERIAL(STRESS,MSTATEV,DDSDDE,STRAN, DSTRAN,
     1 TIME,DTIME,NTENS,NSTATV,PROPS,NPROPS,COORDS)

      INCLUDE 'aba_param.inc'

      CHARACTER*80 CMNAME

C
C      REAL :: equivStress, fractureWorkInc, smean, stressPower, two
      REAL(KIND=8), DIMENSION(1:NTENS)::STRESS,STRAN,DSTRAN
      REAL(KIND=8), DIMENSION(1:NSTATV)::MSTATEV
      REAL(KIND=8), DIMENSION(NTENS,NTENS)::DDSDDE
      REAL(KIND=8), DIMENSION(1:2)::TIME(2),COORDS

      REAL(KIND=8), PARAMETER :: PI = 3.1415926535897932384626433832795d0
      INTEGER, PARAMETER ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      COMMON /KCOMMON_Vars/young, poisson,k_1,th_del,unit_conv
      REAL(KIND=8) :: young, poisson, k_1, th_del, unit_conv
      REAL(KIND=8)  :: zeta, sum_bulk, bulk, zeta0
      COMMON /KM7f_Microplane_1/ qn, ql, qm, w 
      REAL(KIND=8), DIMENSION(1:6,1:np) :: qn, ql, qm
      REAL(KIND=8), DIMENSION(1:np) :: w
      REAL(KIND=8) , DIMENSION(1:6,1:np) :: E_tan_vol, E_tan_dev, E_tan_norm, 
     $    E_tan_shear, E_tan_ela, E_tan_N, E_tan_L, E_tan_M
      REAL(KIND=8) , DIMENSION(1:np) :: deps_N, deps_L, deps_M, 
     $     eps_N, eps_L, eps_M,  sig_N
      REAL(KIND=8) ,DIMENSION(1:np) :: sN_ini,ded,ed_ini,ed_fin
      REAL(KIND=8) , DIMENSION(1:np) :: eN_fin, sN_ela, sN_fin,
     $     sN_boundary, sd_fin, sdneg, sdpos
      REAL(KIND=8) , DIMENSION(1:np) :: sL_ini, sM_ini, sL_fin,
     $     sM_fin, deL, deM, eps_N0_neg, eps_N0_pos

      REAL(KIND=8) , DIMENSION(1:nvhi+2) :: vh_ini, vh_fin 
      REAL(KIND=8) , DIMENSION(1:np) :: E_N
      SAVE    :: /KCOMMON_Vars/
      SAVE    :: /KM7f_Microplane_1/
      
      REAL(KIND=8),SAVE :: k_1i

      
      REAL(KIND=8) , DIMENSION(1:6) :: eps, deps, sig_old 
      REAL(KIND=8) , DIMENSION(1:6) :: sig
      REAL(KIND=8) , DIMENSION(1:6,1:6) :: jacobian
      REAL(KIND=8)  :: sv_ini, ev_ini, dev, ev_fin, sum_sN_fin, devv
      REAL(KIND=8)  :: sv_fin, phi0, phi
      REAL(KIND=8)  :: svneg, ran_num
      INTEGER  ::  i, IC, marks(1:np), rs_size
      
C     ESTABLISH THE SYSTEM
      IF (TIME(2) <= DTIME) THEN
          CALL inputparams()
          CALL setsystem()
C		  CALL random_number(ran_num)
c          k_1i = k_1!(0.95d0+0.025d0*ran_num)*k_1
c          write(6,*), "Mk_1_704=",k_1
          MSTATEV(190) = k_1
      END IF

c      k_1 = MSTATEV(190)
c      write(6,*), "Mk_1_709=",k_1
c      write(6,*), "MSTATEV=",MSTATEV
      vh_ini = MSTATEV 
     
      IF (TIME(2) <= DTIME) THEN
          vh_ini(2) = 1.d0
      END IF
C     Note STRAN and DSTRAN in this code were written in engineering notation
C     Shear strains are shear angles in microplane model
c      write(*,*) 'sig_old is ', kind(sig_old)
c      write(*,*) 'STRESS is ', kind(STRESS)
      sig_old =0.0d0

      eps = 0.0d0
      deps = 0.0d0
      sig_old = [STRESS(1), STRESS(2), 0.D0,STRESS(3), 0.D0,0.D0]*unit_conv
c      write(6,*), "STRESS_OLD=",STRESS
c      write(6,*), "MSTRAN=",STRAN
      eps = [STRAN(1), STRAN(2), 0.D0 ,STRAN(3), 0.D0, 0.D0]
c      write(6,*), "Meps=",eps
      deps = [DSTRAN(1), DSTRAN(2), 0.D0, DSTRAN(3), 0.D0, 0.D0]
c      write(6,*), "Mdeps=",deps
      

             
C     -------------------------------------------------------------------
C     Compute material response :
C     Bounding curve formulation for volumetric,
C ... deviatoric and shear variables.
C     Shear strains are shear angles.
C     -------------------------------------------------------------------

C     ---------------
C ... Initializations
C     ---------------
      sig = 0.d0
      E_tan_N = 0.d0
      E_tan_L = 0.d0
      E_tan_M = 0.d0
C     --------------------------------------------
C ... Evaluate local strains (eps_N, eps_L, eps_M)
C     --------------------------------------------
  
       eps_N =  matmul( eps,qn) 
       eps_L =  matmul( eps,ql)
       eps_M =  matmul( eps,qm)
      deps_N =  matmul(deps,qn)
      deps_L =  matmul(deps,ql)
      deps_M =  matmul(deps,qm)
c        write(6,*), "Meps_N=",eps_N
c        write(6,*), "Meps_L=",eps_L
c        write(6,*), "Meps_M=",eps_M
C      IF (maxval(eps_N) > th_del) THEN
C         vh_fin(188) = 0.d0
C      ELSE
         vh_fin(188) = 1.d0
C      END IF
C----------------------------------------------
C     Read STATEV and STRAN
C----------------------------------------------       
      sv_ini = vh_ini(1)
      phi0 = vh_ini(2)
      sN_ini = vh_ini(3:nvhi:nvhm)
      eps_N0_pos = vh_ini(6:nvhi:nvhm)
      eps_N0_neg = vh_ini(7:nvhi:nvhm)

      zeta0 = vh_ini(189)
      zeta = zeta0      
      
      ev_ini = (eps(1)+ eps(2)+ eps(3))/3.d0
      dev = (deps(1)+deps(2)+deps(3))/3.d0
C----------------------------------------------
C     Normal ELASTIC prediction
C----------------------------------------------      
      CALL c_norm_elastic(eps_N, deps_N, sN_ini, eps_N0_pos,
     $     eps_N0_neg, sv_ini, sN_ela, E_N, zeta, E_tan_ela)

C----------------------------------------------
C ... Normal volumetric law (same for all microplanes)
C----------------------------------------------

      ev_fin = ev_ini + dev 

      CALL c_vol(dev, ev_ini, sv_ini, deps_N, eps_N,svneg,E_tan_vol)

C----------------------------------------------
C ... Normal deviatoric law
C----------------------------------------------
      
      ded = deps_N - dev
      ed_ini = eps_N - ev_ini     
      ed_fin = ed_ini + ded

      CALL c_dev(ded, ed_ini, dev, ev_ini, sv_ini, sd_fin, 
     $     sdneg, sdpos, E_tan_dev)

      eN_fin = ev_fin + ed_fin 

C----------------------------------------------
C ... Normal associated law
C----------------------------------------------      
      CALL c_norm(eN_fin, sv_ini, sN_boundary, E_tan_norm)
      phi0 = vh_ini(2)

      DO i = 1, np
          IF (sN_ela(i) > sN_boundary(i)) THEN 
              marks(i) = 1
              sig_N(i) = sN_boundary(i)
              eps_N0_pos(i) = eN_fin(i)
              E_tan_N(1:6,i) = E_tan_norm(1:6,i)
          ELSE IF (sN_ela(i) < svneg+sdneg(i)) THEN 
              marks(i) = 2
              sig_N(i) = svneg+sdneg(i)
              eps_N0_neg(i) = eN_fin(i)
              E_tan_N(1:6,i) = E_tan_vol(1:6,i) + E_tan_dev(1:6,i)
          ELSE
              marks(i) = 0
              sig_N(i) = sN_ela(i)
              E_tan_N(1:6,i) = E_tan_ela(1:6,i)
          END IF
      END DO

      sum_sN_fin = dot_product(sig_N,w)              
      sv_fin = sum_sN_fin/3.d0 

      sum_bulk = dot_product(E_N,w)   
      bulk = sum_bulk/3.d0

      IF(sv_ini > 0.d0) THEN
          IF(sv_fin > 0.d0) THEN
              devv = abs(dev - ((sv_fin-sv_ini)/bulk))
              zeta = zeta0 + abs(devv)
          END IF
      END IF
C----------------------------------------------
C ... History variables
C----------------------------------------------
      sL_ini = vh_ini(4:nvhi:nvhm) 
      sM_ini = vh_ini(5:nvhi:nvhm) 
 
      sN_fin = sig_N

      eN_fin = eps_N + deps_N
      deL = deps_L
      deM = deps_M

C----------------------------------------------
C ... Shear law
C----------------------------------------------

      CALL c_shear2(eps_L, eps_M, sN_fin, deL, deM, sL_ini
     $     , sM_ini, E_tan_N, sL_fin, sM_fin, ev_fin, E_tan_L, E_tan_M)
c      write(6,*), "sL_ini"
c      write(6,*), sL_ini
c      write(6,*),"00000000000000000000000000000"     
c      write(6,*), "sM_ini"
c      write(6,*), sM_ini
c      write(6,*),"00000000000000000000000000000"          
c      write(6,*), "E_tan_M"
c      write(6,*), E_tan_M
c      write(6,*),"00000000000000000000000000000"
C----------------------------------------------
C     Final Stress Vector
C----------------------------------------------
      sig = matmul(qn,sN_fin*w) + matmul(qm,sM_fin*w) + 
     $        matmul(ql,sL_fin*w)
      DO IC = 1, 6
          E_tan_N(IC,1:np) = E_tan_N(IC,1:np)*w
          E_tan_L(IC,1:np) = E_tan_L(IC,1:np)*w
          E_tan_M(IC,1:np) = E_tan_M(IC,1:np)*w
      END DO
c      write(6,*),"E_tan_N", E_tan_N
c      write(6,*),"E_tan_M", E_tan_M
c      write(6,*),"E_tan_L", E_tan_L
      jacobian = matmul(qn,TRANSPOSE(E_tan_N)) + 
     $ matmul(qm,TRANSPOSE(E_tan_M)) + matmul(ql,TRANSPOSE(E_tan_L))

      sig = sig/unit_conv
      jacobian = jacobian/unit_conv
C
C----------------------------------------------
C ... Update microplane normal and shear stresses
C----------------------------------------------
      
      vh_fin(1) = sv_fin

      vh_fin(3:nvhi:nvhm) = sN_fin 
      vh_fin(4:nvhi:nvhm) = sL_fin 
      vh_fin(5:nvhi:nvhm) = sM_fin 
      vh_fin(6:nvhi:nvhm) = eps_N0_pos
c      write(6,*), "eps_N0_neg"
c      write(6,*), eps_N0_neg
c      write(6,*),"00000000000000000000000000000"      
      vh_fin(7:nvhi:nvhm) = eps_N0_neg 
      vh_fin(189) = zeta
c      write(6,*), "vh_fin=",vh_fin
      MSTATEV = vh_fin
      MSTATEV(190) = k_1
c      write(6,*), "MSTATEV=",MSTATEV
c      write(6,*), "k_1=",k_1
      STRESS = [sig(1), sig(2), sig(4)]
c      write(6,*),"MSTRESS", STRESS
      DDSDDE = RESHAPE(
     $[jacobian(1,1),jacobian(1,2),jacobian(1,4),
     $jacobian(2,1),jacobian(2,2),jacobian(2,4),
     $jacobian(4,1),jacobian(4,2),jacobian(4,4)]
     $,[3,3])
c      IF (NOEL == 454) THEN
c      write(6,*), "MTIME=",TIME(2)
c      write(6,*), "Mjacobian=",jacobian
c      write(6,*), "MDDSDDE=",DDSDDE
c      write(6,*), "Meps=",eps
c      write(6,*), "MCOORDS=",COORDS
c      write(6,*), "Msig_old=",sig_old    
c      write(6,*), "MNTENS=",NTENS  
c      END IF  
C      write(6,*) jacobian
      ! STATEV and STRESS does not update every iteration, but at each increment
      ! Therefore, in each iteration, you cannot update them unless converge
      


C      smean = ( STRESS(1) + STRESS(2) + STRESS(3) )/3.D0
C      equivStress = sqrt( 3.D0/2.D0 * ( (STRESS(1)-smean)**2 +
C     $        (STRESS(2)-smean)**2 + (STRESS(3)-smean)**2 +
C     $        2.D0 * STRESS(4)**2 ) )


      RETURN
      END SUBROUTINE M7FMATERIAL     

C +--------------------------------------------------------------------+
C |                 SUBROUTINE C_NORM_ELASTIC                          |
C +--------------------------------------------------------------------+
      SUBROUTINE c_norm_elastic(eN_ini,deN,sN_ini,eps_N0_pos,eps_N0_neg
     $ ,sv_ini, sN_ela, E_N, zeta,c_tan_ela) 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL(KIND=8), PARAMETER :: PI=3.1415926535897932384626433832795d0
      INTEGER, PARAMETER ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      COMMON/KCOMMON_Vars/young, poisson, k_1, th_del,unit_conv
      REAL(KIND=8) :: young, poisson, k_1, th_del, unit_conv
      COMMON /KM7fBounds_M7fIO/ k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8, c_9, c_10, c_11,c_12,
     $     c_13,c_14,c_15,c_16, c_17,c_18,c_19,c_20,c_21
      REAL(KIND=8) :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9
      REAL(KIND=8) :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8, c_9, c_10
      REAL(KIND=8) :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18, c_19,c_20,c_21
      COMMON /KM7f_Microplane_1/ qn, ql, qm, w 
      REAL(KIND=8), DIMENSION(1:6,1:np) :: qn, ql, qm
      REAL(KIND=8), DIMENSION(1:np)     :: w
      SAVE    :: /KM7f_Microplane_1/
      SAVE :: /KCOMMON_Vars/, /KM7fBounds_M7fIO/
      REAL(KIND=8), DIMENSION(1:np) :: eN_ini, deN, sN_ini
      REAL(KIND=8), DIMENSION(1:np) :: eps_N0_pos, eps_N0_neg, sN_ela, E_N
      REAL(KIND=8), DIMENSION(1:6,1:np) :: c_tan_ela
      REAL(KIND=8) :: sv_ini, zeta

      REAL(KIND=8)  :: E_N_0
      REAL(KIND=8)  :: weight_factor,t0,t1,t2,t3,t4,tsum,fzeta
      INTEGER  :: i 
      
      E_N_0 = young / (1.d0 - 2.d0*poisson)


      DO i=1,np
          IF (sN_ini(i) >= 0.d0) THEN 
              weight_factor = 4.0d1
              t0 = 1.d0
              t1 = (weight_factor*zeta)**2.d0
              t2 = (weight_factor*zeta)**4.d0
              t3 = (weight_factor*zeta)**6.d0
              t4 = (weight_factor*zeta)**8.d0
              tsum = t0+t1+t2
              fzeta = 1.d0/tsum
              E_N(i) = E_N_0*fzeta*exp(-c_19*eps_N0_pos(i))
              IF (sN_ini(i) > E_N_0*(eN_ini(i) + deN(i)) .and. 
     $            sN_ini(i)*deN(i) < 0.d0) THEN 
                  E_N(i) = E_N_0                
              END IF                         
          ELSE 
              E_N(i) = E_N_0*(exp(-c_20*abs(eps_N0_neg(i))/
     $           (1.d0+c_18*max(-sv_ini,0.d0)/E_N_0))+
     $           c_21*max(-sv_ini,0.d0)/E_N_0)
          END IF
          c_tan_ela(1:6,i) = E_N(i)*qn(1:6,i)
      END DO
      sN_ela = sN_ini + E_N*deN
      RETURN
      END SUBROUTINE c_norm_elastic      
      
C +--------------------------------------------------------------------+
C |                        SUBROUTINE C_VOL                            |
C +--------------------------------------------------------------------+
      SUBROUTINE c_vol(dev,ev_ini,sv_ini,deps_N,eps_N,svneg,c_tan_vol) 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL(KIND=8), PARAMETER :: PI=3.1415926535897932384626433832795d0
      REAL(KIND=8), PARAMETER :: sv_ini_p=250.d0
      INTEGER, PARAMETER ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      COMMON /KCOMMON_Vars/young, poisson,k_1,th_del,unit_conv  
      REAL(KIND=8) :: young, poisson, k_1, th_del, unit_conv
      COMMON /KM7fBounds_M7fIO/ k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21
      REAL(KIND=8) :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9
      REAL(KIND=8) :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8, c_9, c_10
      REAL(KIND=8) :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,c_19,c_20,c_21
      COMMON /KM7f_Microplane_1/ qn, ql, qm, w 
      REAL(KIND=8), DIMENSION(1:6,1:np) :: qn, ql, qm
      REAL(KIND=8), DIMENSION(1:np)     :: w
      SAVE    :: /KM7f_Microplane_1/
      SAVE :: /KCOMMON_Vars/, /KM7fBounds_M7fIO/
      
      REAL(KIND=8):: dev,ev_ini,sv_ini,svneg
      REAL(KIND=8)  :: Cv0,ev_fin,xk0, prStrainDIFf,xk4,beta_e     
      REAL(KIND=8), DIMENSION(1:np) :: eps_N, deps_N
      REAL(KIND=8), DIMENSION(1:6,1:np) :: c_tan_vol
      REAL(KIND=8) , DIMENSION(1:6) :: unit_vector
      INTEGER  :: i
      unit_vector = (/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)
      xk0  = k_3*k_1*young 
      Cv0  = young / (1.d0 - 2.d0 * poisson) 
      prStrainDIFf = maxval(eps_N+deps_N,1) - minval(eps_N+deps_N,1)
      
C      write(6,*) k_1
      
      xk4 = (k_6*(prStrainDIFf/k_1)**k_7)
     $     /(1.d0 + min(max(-sv_ini,0.d0),sv_ini_p)/Cv0) + k_4
      
      ev_fin = ev_ini + dev 
      svneg = -xk0*exp(-ev_fin/(xk4*k_1)) 
      
      DO i=1,np
          IF (eps_N(i)+deps_N(i) == maxval(eps_N+deps_N,1)) THEN
              beta_e = 1.d0
          ELSE IF (eps_N(i)+deps_N(i) == minval(eps_N+deps_N,1)) THEN
              beta_e = -1.d0
          ELSE
              beta_e = 0.d0
          END IF
          c_tan_vol(1:6,i) = xk0*exp(-ev_fin/(xk4*k_1))*
     $    (1.d0/(3.d0*xk4*k_1)*unit_vector - ev_fin/k_1*1.d0/xk4**2.d0*
     $    k_6/(1.d0 + min(max(-sv_ini,0.d0),sv_ini_p)/Cv0)*
     $    k_7*(prStrainDIFf/k_1)**(k_7-1.d0)*beta_e/k_1*qn(1:6,i))
      END DO

      END SUBROUTINE c_vol
       
! +--------------------------------------------------------------------+
! |                      SUBROUTINE C_DEV                              |
! +--------------------------------------------------------------------+
      SUBROUTINE c_dev(ded, ed_ini, dev, ev_ini, sv_ini,  sd_fin,
     $     sdneg, sdpos, c_tan_dev) 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL(KIND=8), PARAMETER :: PI=3.1415926535897932384626433832795d0
      INTEGER, PARAMETER ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      COMMON /KCOMMON_Vars/young,poisson,k_1, th_del,unit_conv
      REAL(KIND=8) :: young, poisson, k_1, th_del, unit_conv
      COMMON /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21
      REAL(KIND=8) :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9
      REAL(KIND=8) :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      REAL(KIND=8) :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21
      COMMON /KM7f_Microplane_1/ qn, ql, qm, w 
      REAL(KIND=8), DIMENSION(1:6,1:np) :: qn, ql, qm
      REAL(KIND=8), DIMENSION(1:np)     :: w
      SAVE    :: /KM7f_Microplane_1/
      SAVE :: /KCOMMON_Vars/, /KM7fBounds_M7fIO/
      REAL(KIND=8), DIMENSION(1:np) :: ded, ed_ini
      REAL(KIND=8), DIMENSION(1:np) :: sd_fin
      REAL(KIND=8), DIMENSION(1:np) :: sdneg, sdpos
      REAL(KIND=8), DIMENSION(1:6,1:np) :: c_tan_dev
      REAL(KIND=8) :: dev, ev_ini, sv_ini

      REAL(KIND=8)  :: c_5_0,c_5_1,c_5_M4,c_6_0,c_6_1,c_6_2,c_6_M4
      REAL(KIND=8)  :: c_7_0,c_7_1,c_7_M4,c_8_0,c_8_1,c_8_M4
      REAL(KIND=8)  :: c_9_0,c_9_1,c_9_2,c_9_M4
      REAL(KIND=8)  :: ev_fin, f_c0, E_0, f_cp, c_40, beta_15, 
     $    beta_19, beta_16, beta_17, beta_18
      REAL(KIND=8) , DIMENSION(1:np) :: ed_fin
      REAL(KIND=8) , DIMENSION(1:6) :: unit_vector
      INTEGER  :: i
c                write(6,*) "---k_1_c_dev---"
c                write(6,*) k_1
c                write(6,*) "-------------"   
c                write(6,*) "---young_c_dev---"
c                write(6,*) young
c                write(6,*) "-------------"       
      unit_vector = (/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)
      ev_fin = ev_ini + dev

      c_5_0 = 1.3d-2
      c_5_1 = 4.d0 
      c_7_0  = 1.2d-2 
      c_7_1 = 35.0d2 
      c_8_0  = 1.2d-2 
      c_8_1 = 20.d0 
      c_6_0 = 4.0d2 
      c_6_1  = 4.0d1 
      c_6_2  = 13.d0
      c_9_0 = 4.0d2 
      c_9_1  = 4.0d1 
      c_9_2  = 13.d0
      c_5_M4 = 3.d0    
      c_6_M4 = 1.30d0 
      c_7_M4 = 10.d-1   
      c_8_M4 = 8.d0  
      c_9_M4 = 0.d-1       
      
      f_c0 = 15.08d0
      E_0 = 20000.d0
      c_40 = 1.0d+0

      f_cp = 90.3d0
      beta_15 = c_5_1*exp(-c_40*(f_cp/young-f_c0/E_0))
      beta_16 = c_8_1*exp(-c_40*(f_cp/young-f_c0/E_0))
      beta_17 = c_7_1*exp(-c_40*(f_cp/young-f_c0/E_0))
      beta_18 = c_6_0*exp(-c_40*(f_cp/young-f_c0/E_0))
      beta_19 = c_9_0*exp(-c_40*(f_cp/young-f_c0/E_0))

      c_5 = beta_15*tanh(c_5_0*max(-ev_fin,0.d0)/k_1) + c_5_M4

      c_8 = beta_16*tanh(c_8_0*max(-ev_fin,0.d0)/k_1) + c_8_M4

      c_7 = beta_17*tanh(c_7_0*max(-ev_fin,0.d0)/k_1) + c_7_M4 

      IF (beta_18*max(-ev_fin/k_1-c_6_1,0.d0)>=log(c_6_2)) THEN
          c_6 = c_6_M4*c_6_2
      ELSE
          c_6 = c_6_M4*exp(beta_18*max(-ev_fin/k_1-c_6_1,0.d0))
      END IF
      
      IF (beta_19*max(-ev_fin/k_1-c_9_1,0.d0)>=log(c_9_2)) THEN
          c_9 = c_9_M4*c_9_2
          beta_e = 0.d0
      ELSE
          c_9 = c_9_M4*exp(beta_19*max(-ev_fin/k_1-c_9_1,0.d0))
          beta_e = c_9_M4*exp(beta_19*max(-ev_fin/k_1-c_9_1,0.d0))*
     $ beta_19*hside(-ev_fin/k_1-c_9_1)*(-1.d0/(3.d0*k_1))  
      END IF

      ed_fin = ed_ini + ded 

      sdneg = -young*k_1*c_8/
     $    (1.d0+(max(-ed_fin-c_9*c_8*k_1,0.d0)/(c_7*k_1))**2.d0) 
      sdpos = young*k_1*c_5/
     $    (1.d0+(max(ed_fin-c_6*c_5*k_1,0.d0)/(c_7*k_1*c_20))**2.d0)
      
      DO i = 1,np
          c_tan_dev(1:6,i) = young*k_1*c_8/
     $ (1.d0+(max(-ed_fin(i)-c_9*c_8*k_1,0.d0)/(c_7*k_1))**2.d0)**2.d0*
     $ 2.d0*max(-ed_fin(i)-c_9*c_8*k_1,0.d0)*
     $ hside(-ed_fin(i)-c_9*c_8*k_1)/(c_7*k_1)**2.d0*
     $ (-qn(1:6,i) + unit_vector/3.d0 - k_1*c_9*beta_16*
     $ (1.d0-(tanh(c_8_0*max(-ev_fin,0.d0)/k_1))**2.d0)*c_8_0/k_1*
     $ hside(-ev_fin)*(-unit_vector/3.d0)-k_1*c_8*beta_e*unit_vector) - 
     $ young*k_1/(1.d0+(max(-ed_fin(i)-c_9*c_8*k_1,0.d0)/(c_7*k_1))**
     $ 2.d0)*beta_16*(1.d0-(tanh(c_8_0*max(-ev_fin,0.d0)/k_1))**2.d0)*
     $ c_8_0/k_1*hside(-ev_fin)*(-unit_vector/3.d0) + young*k_1*c_8/
     $ (1.d0+(max(-ed_fin(i)-c_9*c_8*k_1,0.d0)/(c_7*k_1))**2.d0)**2.d0*
     $ (max(-ed_fin(i)-c_9*c_8*k_1,0.d0)/k_1)**2.d0*(-2.d0/c_7**3.d0)*
     $ beta_17*(1.d0-(tanh(c_7_0*max(-ev_fin,0.d0)/k_1))**2.d0)*
     $ c_7_0/k_1*hside(-ev_fin)*(-unit_vector/3.d0)
      END DO
      
      RETURN
      END SUBROUTINE c_dev  


C +--------------------------------------------------------------------+
C |                         SUBROUTINE C_NORM                          |
C +--------------------------------------------------------------------+
      SUBROUTINE c_norm(eN_fin,sv_ini,sN_boundary,c_tan_norm) 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL(KIND=8), PARAMETER :: PI=3.1415926535897932384626433832795d0
      INTEGER, PARAMETER ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      COMMON /KCOMMON_Vars/young,poisson,k_1, th_del,unit_conv
      REAL(KIND=8) :: young, poisson, k_1, th_del, unit_conv
      COMMON /KM7fBounds_M7fIO/ k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21
      REAL(KIND=8) :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9
      REAL(KIND=8) :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      REAL(KIND=8) :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21
      COMMON /KM7f_Microplane_1/ qn, ql, qm, w 
      REAL(KIND=8), DIMENSION(1:6,1:np) :: qn, ql, qm
      REAL(KIND=8), DIMENSION(1:np)     :: w
      SAVE    :: /KM7f_Microplane_1/
      SAVE :: /KCOMMON_Vars/, /KM7fBounds_M7fIO/
      REAL(KIND=8), DIMENSION(1:np) :: eN_fin, sN_boundary
      REAL(KIND=8), DIMENSION(1:6,1:np) :: c_tan_norm
      REAL(KIND=8) :: sv_ini
      REAL(KIND=8)  :: fstar, beta_N, E_N_0, eb_N
      REAL(KIND=8)  :: d_1,d_2,d_3,d_4,d_5,d_6,V_f
      INTEGER  :: i
c                write(6,*) "---k_1_c_norm---"
c                write(6,*) k_1
c                write(6,*) "-------------"  
      d_1 = 0.095d0 
      d_2 = 35.d0   
      d_3 = 1.7d0   
      d_4 = 1.7d0 
      d_5 = 1.d3
      d_6 = 25.d0
      V_f = 0.d0
      E_N_0 = young/(1.d0 - 2.d0*poisson)
      c_1 = d_1*tanh(d_2*V_f-d_3) + 
     $    d_4*exp(-max(-sv_ini-d_6,0.d0)/E_N_0*d_5) 

      IF (sv_ini .lt. 0.d0) THEN 
         eb_N = c_3*k_1 - c_4 / E_N_0 * sv_ini  
      ELSE 
         eb_N = c_3*k_1
      END IF       

      fstar  = k_1*young*c_1 
      beta_N = c_2*c_1*k_1 

      sN_boundary = fstar*exp(-max(eN_fin-beta_N,0.d0)/eb_N)
      
      DO i = 1,np
      c_tan_norm(1:6,i) = fstar*exp(-max(eN_fin(i)-beta_N,0.d0)/eb_N)*
     $ (-hside(eN_fin(i)-beta_N)/eb_N)*qn(1:6,i)
      END DO

      
      RETURN 
      END SUBROUTINE c_norm



! +--------------------------------------------------------------------+
! |                          SUBROUTINE C_SHEAR2                       |
! +--------------------------------------------------------------------+
      SUBROUTINE c_shear2(eps_L, eps_M, sN_fin, deL, deM,
     $      sL_ini,sM_ini,c_tan_N,sL_fin,sM_fin,ev_fin,c_tan_L,c_tan_M)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL(KIND=8), PARAMETER :: PI=3.1415926535897932384626433832795d0
      INTEGER, PARAMETER ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      COMMON /KCOMMON_Vars/ young,poisson,k_1, th_del,unit_conv
      REAL(KIND=8) :: young, poisson, k_1, th_del, unit_conv
      COMMON /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21
      REAL(KIND=8) :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9
      REAL(KIND=8) :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,c_9, c_10
      REAL(KIND=8) :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21
      COMMON /KM7f_Microplane_1/ qn, ql, qm, w 
      REAL(KIND=8), DIMENSION(1:6,1:np) :: qn, ql, qm
      REAL(KIND=8), DIMENSION(1:np)     :: w
      SAVE    :: /KM7f_Microplane_1/
      SAVE :: /KCOMMON_Vars/, /KM7fBounds_M7fIO/
      REAL(KIND=8), DIMENSION(1:np)  :: eps_L, eps_M, sN_fin, deL,
     $     deM, sL_ini, sM_ini, sL_fin, sM_fin 
      REAL(KIND=8), DIMENSION(1:6,1:np) :: c_tan_N, c_tan_L, c_tan_M
      REAL(KIND=8) :: ev_fin
      REAL(KIND=8)  :: c_12_p, c_13_p, c_6_p, s_0, sig_o, Ct, sT
      REAL(KIND=8) , DIMENSION(1:np) :: E_T, fsp, stfe, fsp_0,
     $ sL_fin_e, sM_fin_e
      REAL(KIND=8) , DIMENSION(1:6) :: unit_vector
      INTEGER  :: i
      fsp = 0.d0
      unit_vector = (/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)
      sT = 0.1d0*young*k_1
      
      
      Ct = young/(1.d0+poisson)*(1.d0-4.d0*poisson)/(1.d0-2.d0*poisson)
      c_6_p  = c_10
      c_12_p = c_11 
      c_13_p = c_12*3.371d-4
      s_0 = k_1*k_2*Ct
c                write(6,*) "---k_1---"
c                write(6,*) k_1
c                write(6,*) "-------------"        
c                write(6,*) "---k_2---"
c                write(6,*) k_2
c                write(6,*) "-------------"                   
c                write(6,*) "---Ct---"
c                write(6,*) Ct
c                write(6,*) "-------------"  
      sig_o = max(Ct*k_1*(c_12_p-c_13_p*max(ev_fin,0.d0)/k_1),0.d0)
c                write(6,*) "---0.5d0*(1-tanh(sN_fin/sT))---"
c                write(6,*) 0.5d0*(1-tanh(sN_fin/sT))
c                write(6,*) "-------------"  
c                write(6,*) "---((c_6_p*max(-sN_fin+sig_o,0.d0))---"
c                write(6,*) (c_6_p*max(-sN_fin+sig_o,0.d0))
c                write(6,*) "-------------" 
c                write(6,*) "---(-1.d0)+s_0**(-1.d0)---"
c                write(6,*) (-1.d0)+s_0**(-1.d0)
c                write(6,*) "-------------"                 
c                write(6,*) "---(1+tanh(sN_fin/sT))---"
c                write(6,*) (1+tanh(sN_fin/sT))
c                write(6,*) "-------------"  
c                write(6,*) "-(c_6_p*max(sig_o, 0.d0))**(-1.d0)+s_0**(-1.d0)--"
c                write(6,*) (c_6_p*max(sig_o, 0.d0))**(-1.d0)+s_0**(-1.d0)
c                write(6,*) "-------------"                                  
      fsp = 0.5d0*(1-tanh(sN_fin/sT))*((c_6_p*max(-sN_fin+sig_o,0.d0))**
     $    (-1.d0)+s_0**(-1.d0))**(-1.d0) + 0.5d0*(1+tanh(sN_fin/sT))*
     $    ((c_6_p*max(sig_o, 0.d0))**(-1.d0)+s_0**(-1.d0))**(-1.d0)
c                write(6,*) "---young---"
c                write(6,*) young
c                write(6,*) "-------------"        
c                write(6,*) "---E_T----"
c                write(6,*) E_T
c                write(6,*) "-------------"     
c                write(6,*) "---sL_ini----"
c                write(6,*) sL_ini
c                write(6,*) "-------------"   
c                write(6,*) "---deL----"
c                write(6,*) deL
c                write(6,*) "-------------"  
c                write(6,*) "---sM_ini----"
c                write(6,*) sM_ini
c                write(6,*) "-------------"    
c                write(6,*) "---deM----"
c                write(6,*) deM
c                write(6,*) "-------------"    
      E_T = Ct
      sL_fin_e = sL_ini + E_T*deL
      sM_fin_e = sM_ini + E_T*deM
      stfe = sqrt(sL_fin_e*sL_fin_e + sM_fin_e*sM_fin_e)
      DO i = 1, np
c                    write(6,*) "---(i)----"
c                    write(6,*) i
c                    write(6,*) "-------------"  
c                    write(6,*) "---stfe(i)----"
c                    write(6,*) stfe(i)
c                    write(6,*) "-------------"   
c                    write(6,*) "---fsp(i) ----"
c                    write(6,*) fsp(i) 
c                    write(6,*) "-------------" 
c                    write(6,*) "---sig_o  ----"
c                    write(6,*) sig_o 
c                    write(6,*) "-------------"  
c                    write(6,*) "---sN_fin(i)  ----"
c                    write(6,*) sN_fin(i)
c                    write(6,*) "-------------"   
          IF (stfe(i) > fsp(i) .and. stfe(i) /= 0.d0) THEN
                  
              IF (sig_o > sN_fin(i) .and. sig_o > 0.d0) THEN
                  sL_fin(i) = sL_fin_e(i)/stfe(i)*fsp(i)
                  sM_fin(i) = sM_fin_e(i)/stfe(i)*fsp(i)
                  c_tan_L(1:6,i)=fsp(i)/stfe(i)*E_T(i)*ql(1:6,i)-fsp(i)*
     $ sL_fin_e(i)/stfe(i)**2.d0*(sL_fin_e(i)/stfe(i)*E_T(i)*ql(1:6,i) +
     $ sM_fin_e(i)/stfe(i)*E_T(i)*qm(1:6,i)) + sL_fin_e(i)/stfe(i)*( 
     $ -0.5d0*(1.d0-(tanh(sN_fin(i)/sT))**2.d0)*(((c_6_p*max(-sN_fin(i)+
     $ sig_o,0.d0))**(-1.d0)+s_0**(-1.d0))**(-1.d0))*c_tan_N(1:6,i)/sT +
     $ 0.5d0*(1.d0-tanh(sN_fin(i)/sT))*(((c_6_p*max(-sN_fin(i)+sig_o,
     $ 0.d0))**(-1.d0)+s_0**(-1.d0))**(-2.d0)*(c_6_p*max(-sN_fin(i)+
     $ sig_o,0.d0))**(-2.d0)*c_6_p*hside(-sN_fin(i)+sig_o)*(hside(Ct*
     $ k_1*(c_12_p-c_13_p*max(ev_fin,0.d0)/k_1))*(-Ct*c_13_p*
     $ hside(ev_fin))*unit_vector/3.d0-c_tan_N(1:6,i))) + 0.5d0*
     $ (1.d0+(tanh(sN_fin(i)/sT))**2.d0)*(((c_6_p*max(sig_o,0.d0))**
     $ (-1.d0)+s_0**(-1.d0))**(-1.d0))*c_tan_N(1:6,i)/sT + 0.5d0*(1.d0+
     $ tanh(sN_fin(i)/sT))*(((c_6_p*max(sig_o,0.d0))**(-1.d0)+s_0**
     $ (-1.d0))**(-2.d0)*(c_6_p*max(sig_o,0.d0))**(-2.d0)*c_6_p*
     $ hside(sig_o)*(hside(Ct*k_1*(c_12_p-c_13_p*max(ev_fin,0.d0)/k_1))*
     $ (-Ct*c_13_p*hside(ev_fin))*unit_vector/3.d0)) )
                       
                  c_tan_M(1:6,i)=fsp(i)/stfe(i)*E_T(i)*qm(1:6,i)-fsp(i)*
     $ sM_fin_e(i)/stfe(i)**2.d0*(sM_fin_e(i)/stfe(i)*E_T(i)*qm(1:6,i) +
     $ sL_fin_e(i)/stfe(i)*E_T(i)*ql(1:6,i)) + sM_fin_e(i)/stfe(i)*( 
     $ -0.5d0*(1.d0-(tanh(sN_fin(i)/sT))**2.d0)*(((c_6_p*max(-sN_fin(i)+
     $ sig_o,0.d0))**(-1.d0)+s_0**(-1.d0))**(-1.d0))*c_tan_N(1:6,i)/sT +
     $ 0.5d0*(1.d0-tanh(sN_fin(i)/sT))*(((c_6_p*max(-sN_fin(i)+sig_o,
     $ 0.d0))**(-1.d0)+s_0**(-1.d0))**(-2.d0)*(c_6_p*max(-sN_fin(i)+
     $ sig_o,0.d0))**(-2.d0)*c_6_p*hside(-sN_fin(i)+sig_o)*(hside(Ct*
     $ k_1*(c_12_p-c_13_p*max(ev_fin,0.d0)/k_1))*(-Ct*c_13_p*
     $ hside(ev_fin))*unit_vector/3.d0-c_tan_N(1:6,i))) + 0.5d0*
     $ (1.d0+(tanh(sN_fin(i)/sT))**2.d0)*(((c_6_p*max(sig_o,0.d0))**
     $ (-1.d0)+s_0**(-1.d0))**(-1.d0))*c_tan_N(1:6,i)/sT + 0.5d0*(1.d0+
     $ tanh(sN_fin(i)/sT))*(((c_6_p*max(sig_o,0.d0))**(-1.d0)+s_0**
     $ (-1.d0))**(-2.d0)*(c_6_p*max(sig_o,0.d0))**(-2.d0)*c_6_p*
     $ hside(sig_o)*(hside(Ct*k_1*(c_12_p-c_13_p*max(ev_fin,0.d0)/k_1))*
     $ (-Ct*c_13_p*hside(ev_fin))*unit_vector/3.d0)) )
c                    write(6,*) "---c_tan_M--1------------"
c                    write(6,*), "i" ,i
c                    write(6,*) "---c_tan_M(1:6,i)----"
c                    write(6,*) c_tan_M(1:6,i)
c                    write(6,*) "-------------"                  
c                    write(6,*) "--********---"
c                    write(6,*) qm
c                    write(6,*) "--******---"
              ELSE IF (sig_o > sN_fin(i) .and. sig_o <= 0.d0) THEN
                  sL_fin(i) = sL_fin_e(i)/stfe(i)*fsp(i)
                  sM_fin(i) = sM_fin_e(i)/stfe(i)*fsp(i)
                  c_tan_L(1:6,i)=fsp(i)/stfe(i)*E_T(i)*ql(1:6,i)-fsp(i)*
     $ sL_fin_e(i)/stfe(i)**2.d0*(sL_fin_e(i)/stfe(i)*E_T(i)*ql(1:6,i) +
     $ sM_fin_e(i)/stfe(i)*E_T(i)*qm(1:6,i)) + sL_fin_e(i)/stfe(i)*( 
     $ -0.5d0*(1.d0-(tanh(sN_fin(i)/sT))**2.d0)*(((c_6_p*max(-sN_fin(i)+
     $ sig_o,0.d0))**(-1.d0)+s_0**(-1.d0))**(-1.d0))*c_tan_N(1:6,i)/sT +
     $ 0.5d0*(1.d0-tanh(sN_fin(i)/sT))*(((c_6_p*max(-sN_fin(i)+sig_o,
     $ 0.d0))**(-1.d0)+s_0**(-1.d0))**(-2.d0)*(c_6_p*max(-sN_fin(i)+
     $ sig_o,0.d0))**(-2.d0)*c_6_p*hside(-sN_fin(i)+sig_o)*(hside(Ct*
     $ k_1*(c_12_p-c_13_p*max(ev_fin,0.d0)/k_1))*(-Ct*c_13_p*
     $ hside(ev_fin))*unit_vector/3.d0-c_tan_N(1:6,i))) )
                       
                  c_tan_M(1:6,i)=fsp(i)/stfe(i)*E_T(i)*qm(1:6,i)-fsp(i)*
     $ sM_fin_e(i)/stfe(i)**2.d0*(sM_fin_e(i)/stfe(i)*E_T(i)*qm(1:6,i) +
     $ sL_fin_e(i)/stfe(i)*E_T(i)*ql(1:6,i)) + sM_fin_e(i)/stfe(i)*( 
     $ -0.5d0*(1.d0-(tanh(sN_fin(i)/sT))**2.d0)*(((c_6_p*max(-sN_fin(i)+
     $ sig_o,0.d0))**(-1.d0)+s_0**(-1.d0))**(-1.d0))*c_tan_N(1:6,i)/sT +
     $ 0.5d0*(1.d0-tanh(sN_fin(i)/sT))*(((c_6_p*max(-sN_fin(i)+sig_o,
     $ 0.d0))**(-1.d0)+s_0**(-1.d0))**(-2.d0)*(c_6_p*max(-sN_fin(i)+
     $ sig_o,0.d0))**(-2.d0)*c_6_p*hside(-sN_fin(i)+sig_o)*(hside(Ct*
     $ k_1*(c_12_p-c_13_p*max(ev_fin,0.d0)/k_1))*(-Ct*c_13_p*
     $ hside(ev_fin))*unit_vector/3.d0-c_tan_N(1:6,i))) )
c                write(6,*) "---c_tan_M--2------------"
c                    write(6,*), "i" ,i
c                    write(6,*) "---c_tan_M(1:6,i)----"
c                    write(6,*) c_tan_M(1:6,i)
c                    write(6,*) "--*******-----"                  
c                    write(6,*) "--********---"
c                    write(6,*) qm
c                    write(6,*) "--******---"
              ELSE IF (sig_o <= sN_fin(i) .and. sig_o > 0.d0) THEN  
                  sL_fin(i) = sL_fin_e(i)/stfe(i)*fsp(i)
                  sM_fin(i) = sM_fin_e(i)/stfe(i)*fsp(i)
                  c_tan_L(1:6,i)=fsp(i)/stfe(i)*E_T(i)*ql(1:6,i)-fsp(i)*
     $ sL_fin_e(i)/stfe(i)**2.d0*(sL_fin_e(i)/stfe(i)*E_T(i)*ql(1:6,i) +
     $ sM_fin_e(i)/stfe(i)*E_T(i)*qm(1:6,i)) + sL_fin_e(i)/stfe(i)*
     $ ( 0.5d0*
     $ (1.d0+(tanh(sN_fin(i)/sT))**2.d0)*(((c_6_p*max(sig_o,0.d0))**
     $ (-1.d0)+s_0**(-1.d0))**(-1.d0))*c_tan_N(1:6,i)/sT + 0.5d0*(1.d0+
     $ tanh(sN_fin(i)/sT))*(((c_6_p*max(sig_o,0.d0))**(-1.d0)+s_0**
     $ (-1.d0))**(-2.d0)*(c_6_p*max(sig_o,0.d0))**(-2.d0)*c_6_p*
     $ hside(sig_o)*(hside(Ct*k_1*(c_12_p-c_13_p*max(ev_fin,0.d0)/k_1))*
     $ (-Ct*c_13_p*hside(ev_fin))*unit_vector/3.d0)) )
                       
                  c_tan_M(1:6,i)=fsp(i)/stfe(i)*E_T(i)*qm(1:6,i)-fsp(i)*
     $ sM_fin_e(i)/stfe(i)**2.d0*(sM_fin_e(i)/stfe(i)*E_T(i)*qm(1:6,i) +
     $ sL_fin_e(i)/stfe(i)*E_T(i)*ql(1:6,i)) + sM_fin_e(i)/stfe(i)*
     $ ( 0.5d0*
     $ (1.d0+(tanh(sN_fin(i)/sT))**2.d0)*(((c_6_p*max(sig_o,0.d0))**
     $ (-1.d0)+s_0**(-1.d0))**(-1.d0))*c_tan_N(1:6,i)/sT + 0.5d0*(1.d0+
     $ tanh(sN_fin(i)/sT))*(((c_6_p*max(sig_o,0.d0))**(-1.d0)+s_0**
     $ (-1.d0))**(-2.d0)*(c_6_p*max(sig_o,0.d0))**(-2.d0)*c_6_p*
     $ hside(sig_o)*(hside(Ct*k_1*(c_12_p-c_13_p*max(ev_fin,0.d0)/k_1))*
     $ (-Ct*c_13_p*hside(ev_fin))*unit_vector/3.d0)) )
c                    write(6,*) "---c_tan_M--3------------"
c                    write(6,*), "i" ,i
c                    write(6,*) "---c_tan_M(1:6,i)----"
c                    write(6,*) c_tan_M(1:6,i)
c                    write(6,*) "--00001384-----"                  
c                    write(6,*) "--********---"
c                    write(6,*) qm
c                    write(6,*) "--******---"
              ELSE     
                  sL_fin(i) = sL_fin_e(i)/stfe(i)*fsp(i)
                  sM_fin(i) = sM_fin_e(i)/stfe(i)*fsp(i)
                  c_tan_L(1:6,i)=fsp(i)/stfe(i)*E_T(i)*ql(1:6,i)-fsp(i)*
     $ sL_fin_e(i)/stfe(i)**2.d0*(sL_fin_e(i)/stfe(i)*E_T(i)*ql(1:6,i) +
     $ sM_fin_e(i)/stfe(i)*E_T(i)*qm(1:6,i))
                  c_tan_M(1:6,i)=fsp(i)/stfe(i)*E_T(i)*qm(1:6,i)-fsp(i)*
     $ sM_fin_e(i)/stfe(i)**2.d0*(sM_fin_e(i)/stfe(i)*E_T(i)*qm(1:6,i) +
     $ sL_fin_e(i)/stfe(i)*E_T(i)*ql(1:6,i))
c                    write(6,*) "---c_tan_M--4------------"
c                    write(6,*), "i" ,i
c                    write(6,*) "---c_tan_M(1:6,i)----"
c                    write(6,*) c_tan_M(1:6,i)
c                    write(6,*) "--*********----"                  
c                    write(6,*) "--********---"
c                    write(6,*) qm
c                    write(6,*) "--******---"
              END IF
          ELSE
                  sL_fin(i) = sL_fin_e(i)
                  sM_fin(i) = sM_fin_e(i)
                  c_tan_L(1:6,i) = E_T(i)*ql(1:6,i)
                  c_tan_M(1:6,i) = E_T(i)*qm(1:6,i)
c                  write(6,*) "---c_tan_M--5------------"
c                  write(6,*), "i" ,i
c                  write(6,*) "---c_tan_M(1:6,i)----"
c                  write(6,*) c_tan_M(1:6,i)
c                  write(6,*) "--00001384-----"                  
c                  write(6,*) "---qm-----"
c                  write(6,*) qm
c                  write(6,*) "--00001384-----"
          END IF
      END DO
c      write(*,*) "---sL_fin--------------"
c      write(*,*) sL_fin
c      write(*,*) "---sM_fin--------------"
c      write(*,*) sM_fin
c      write(*,*) "---c_tan_L--------------"
c      write(*,*) c_tan_L
c      write(6,*) "---c_tan_M1390--------------"
c      write(6,*) c_tan_M
c      write(6,*) "------0000000000000-1390-----------"
      
      RETURN 
      END SUBROUTINE c_shear2 

      FUNCTION hside(x) !Heaviside function
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL(KIND=8) :: hside, x
      hside = 1.0D0
      IF (x < 0.0D0 ) hside = 0.0D0
      END FUNCTION hside


C **********************************************************************
C *** SUBROUTINE INPUTPARAMS *******************************************
C **********************************************************************
      SUBROUTINE inputparams() 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL(KIND=8), PARAMETER :: PI=3.1415926535897932384626433832795d0
      INTEGER, PARAMETER :: np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      COMMON /KCOMMON_Vars/young, poisson,k_1, th_del,unit_conv  
      REAL(KIND=8) :: young, poisson, k_1, th_del, unit_conv
      COMMON /KM7fBounds_M7fIO/k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9,
     $     c_1, c_2, c_3, c_4, c_5,
     $     c_6, c_7, c_8, c_9, c_10, c_11,c_12,c_13,c_14,c_15,c_16,
     $     c_17,c_18,c_19,c_20,c_21
      REAL(KIND=8) :: k_2, k_3, k_4, k_5, k_6, k_7, k_8, k_9
      REAL(KIND=8) :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8,
     $     c_9, c_10
      REAL(KIND=8) :: c_11,c_12,c_13,c_14,c_15,c_16,c_17,c_18,
     $     c_19,c_20,c_21
      SAVE :: /KCOMMON_Vars/,/KM7fBounds_M7fIO/

      th_del = 0.005d10 

      unit_conv=1.d0                



      young  = 28700.0d0 
      poisson= 0.18000d0
      k_1    = 85.00d-6 
      k_2    = 110.000d0
      k_3    = 12.d0 
      k_4    = 38.d0 
      

      k_8 = 3.d0
      k_9 = 5.d-1
      k_5 = 1.d0
      k_6 = 1.0d-4 
      k_7 = 1.8d0 
     
      c_2 = 1.76d-1 
      c_3 = 12.d0 
      c_4 = 10.d0 
      c_10 = 3.3d-1 
      c_11 = 5.d-1 
      c_12 = 7.00d3 
      c_16 = 10.d0 
      c_17 = 1.00d-2       
      c_18 = 4.d3 
      c_19 = 4.5d3 
      c_20 = 3.d2 
      c_21 = 6.d1 
     
     
        
      RETURN
      END SUBROUTINE inputparams