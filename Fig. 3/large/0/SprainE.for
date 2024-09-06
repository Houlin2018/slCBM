      subroutine SprainE(Sprainlocal,PROPS,NPROPS,bmata,bmatb,bmatc,nmatb,NDOFEL,
     1 MLVARX,U,DU,D,N,stress,strain,fd,feta,flambda,K_f,COORDS,sp_energy)
      implicit none

      ! Input arrays and variables
      INTEGER ::      NDOFEL,MLVARX
      real(kind=8) :: bmata(6, 24),bmatb(8, 32),bmatc(4, 24),nmatb(4, 32)
      real(kind=8) :: U(NDOFEL),DU(MLVARX,1),D(6, 6),N(8),stress(6),strain(6)
      real(kind=8) :: PROPS(*),COORDS(3,8)
      integer ::       NPROPS
      ! Output arrays
      real(kind=8):: fd(24,1),feta(32,1),flambda(32,1),K_s(88, 88),K_f(88, 88)

      ! Variables
      real(kind=8) :: dd(24,1), eta(32,1), lambda(32,1),nlambda(4,32)
      real(kind=8) :: dd_new(24,1), eta_new(32,1), lambda_new(32,1),dT_new(1,24)
      real(kind=8) :: dT(1,24),etaT(1,32),nlambdaT(32,4)
      real(kind=8) :: kappa, C_p, A(1,1),A0,l_0

      ! Temporary variables
      real(kind=8) :: temp_matrix1(8, 3),temp,mc_result,h_result
      real(kind=8) :: temp_matrix2(4, 1),temp_matrix3(4,1)
      real(kind=8) :: temp_matrix5(1, 1)
!      real(kind=8) :: temp_matrix6(1, 1)
      real(kind=8) :: temp_matrix7(6, 1)
      real(kind=8) :: temp_vector1(8, 1)

      real(kind=8) :: omega2,Dzz(8,8)
      real(kind=8) :: Sprainlocal(2),eqSprain,gamma0,gammac,Spress(8,1)
      real(kind=8) :: sp_energy
      integer(kind=8):: ind_dof1(24),ind_dof2(32),ind_dof3(32)
      integer(kind=8):: i
C      write(6,*), "N shape function=",N	
C      write(6,*), "bmata=",bmata	
C      write(6,*), "bmatb=",bmatb
C      write(6,*), "bmatc=",bmatc
C      write(6,*), "nmatb=",nmatb
C      write(6,*), "U=",U
C      write(6,*), "DU=",DU(:,1)
!      write(6,*), "DU IS",KIND(DU)
!      PAUSE
      fd = 0d0
      feta=0d0
      flambda=0d0
      Dzz = 0.0d0
      ind_dof1 =   (/1,2,3, 12,13,14, 23,24,25, 34,35,36,
     $               45,46,47, 56,57,58, 67,68,69, 78,79,80/)
      ind_dof2 =   (/4,5,6,7,15,16,17,18,26,27,28,29,37,38,39,40,
     $               48,49,50,51,59,60,61,62,70,71,72,73,81,82,83,84/)
      ind_dof3 =   (/8,9,10,11,19,20,21,22,30,31,32,33,41,42,43,44,
     $               52,53,54,55,63,64,65,66,74,75,76,77,85,86,87,88/)
      do i = 1,24
         dd(i,1)=U(ind_dof1(i))
      end do
      dT = transpose(dd)
      do i =1,32
         eta(i,1)=U(ind_dof2(i))
      end do
      etaT = transpose(eta)
      do i =1,32
         lambda(i,1)=U(ind_dof3(i))
      end do

      do i = 1,24
      dd_new(i,1)=U(ind_dof1(i))!+DU(ind_dof1(i),1)
      end do
C      write(6,*), "dd_new=",dd_new
      dT_new = transpose(dd_new)
      do i =1,32
      eta_new(i,1)=U(ind_dof2(i))!+DU(ind_dof2(i),1)
      end do
C      write(6,*), "eta_new=",eta_new
      do i =1,32
      lambda_new(i,1)=U(ind_dof3(i))!+DU(ind_dof3(i),1)
      end do
C      write(6,*), "lambda_new=",lambda_new
      nlambda=0d0
      do i = 1,8
      nlambda(1,1+4*(i-1))=N(i)
      nlambda(2,2+4*(i-1))=N(i)
      nlambda(3,3+4*(i-1))=N(i)
      nlambda(4,4+4*(i-1))=N(i)
      end do
      nlambdaT = transpose(nlambda)


      kappa = 50d0
      l_0   = 5d0

      gamma0 = PROPS(5)
      gammac = PROPS(6)  
C      write(6,*), "gamma0=",gamma0
C      write(6,*), "gammac=",gammac    
      
      



C      write(6,*), "A(1,1)=",A(1,1)
C      write(6,*), "A0 new=",A0
      ! Calculate fd
      temp_matrix7(1:6,1) =stress 
      fd = matmul(transpose(bmata), temp_matrix7)
      fd = fd + matmul(matmul(transpose(bmatc), nlambda),lambda_new)
C      write(6,*), "fd=",fd	
      ! Calculate feta
      call Spress_Sprain(Sprainlocal,eta_new,bmatb,eqSprain,kappa,gamma0,gammac,Spress,Dzz,sp_energy)
      omega2 =Sprainlocal(2)
      sp_energy = l_0**2.0d0*sp_energy
C      write(6,*), "sp_energy=",sp_energy      
      feta = l_0**2.0d0*matmul(transpose(bmatb), Spress)
      feta =feta - matmul(transpose(nmatb),matmul(nlambda,lambda_new))
C      write(6,*), "feta=",feta	
      ! Calculate flambda
      temp_matrix2 = matmul(bmatc,dd_new)
      flambda = matmul(nlambdaT, temp_matrix2)
      temp_matrix3 = matmul(nmatb,eta_new)
      flambda = flambda - matmul( nlambdaT,temp_matrix3)
c      flambda = 0d0
C      write(6,*), "flambda=",flambda
      ! Calculate K_s matrix blocks
      ! K_s matrix blocks:
      ! dd, deta, dlambda
      ! etad, etaeta, etalambda
      ! lambdad, lambdaeta, lambdalambda
c      call Spress_Sprain(Sprainlocal,eta,bmatb,eqSprain,kappa,gamma0,gammac,Spress,Dzz)
      K_s = 0d0
      ! Block dd 
!      temp_matrix6 = matmul(nlambda,lambda)
      K_s(1:24, 1:24) = matmul(matmul(transpose(bmata),D),bmata)!+
c      write(6,*), "bmata in sprianE=",bmata
!     $              2d0* temp_matrix6(1,1)*B3TB3
C      write(6,*), "K_s(1:8, 1:8)=", K_s(1:8, 1:8)
c      write(6,*), "0.125 K_s(1:24, 1:24)=", 0.125d0*K_s(1:24, 1:24)
      ! Block deta 
      K_s(1:24, 25:56) = 0d0
c      write(6,*), "K_s(1:8, 9:24)=", K_s(1:8, 9:24)
      ! Block dlambda 
      K_s(1:24, 57:88) =  matmul(transpose(bmatc),nlambda)
c      write(6,*), "0.125 K_s(1:24, 57:88) dlambda=", 0.125*K_s(1:24, 57:88)
      ! Block etad 
      K_s(25:56, 1:24) = 0d0
c      write(6,*), "K_s(9:24, 1:8)=", K_s(9:24, 1:8)
      ! Block etaeta 
      K_s(25:56, 25:56) = l_0*l_0*matmul(transpose(bmatb),matmul(Dzz,bmatb))
c      write(6,*), "0.125 K_s(25:56, 25:56)=", 0.125*K_s(25:56, 25:56)
      ! Block etalambda 
      K_s(25:56, 57:88) =  -matmul(transpose(nmatb),nlambda)
c      write(6,*), "K_s(25:56, 57:88) etalambda=",0.125*K_s(25:56, 57:88)
      ! Block lambdad  bmatc
      K_s(57:88, 1:24) = matmul(nlambdaT,bmatc)
c      write(6,*), "K_s(57:88, 1:24)=",0.125*K_s(57:88, 1:24)
      ! Block lambdaeta 
      K_s(57:88, 25:56) = - matmul(nlambdaT,nmatb)
c      write(6,*), "K_s(57:88, 25:56) lambdaeta=",0.125*K_s(57:88, 25:56)
      ! Block lambdalambda 
      K_s(57:88, 57:88) = 0.d0
c      write(6,*), "K_s(25:40, 25:40)=",K_s(25:40, 25:40)
c      write(6,*), "K_s=",K_s
      call finalmatrix(K_s,K_f)
c      pause

      end subroutine SprainE

*******************************************************************
      subroutine Spress_Sprain(Sprainlocal,eta,bmatb,eqSprain,kappa0,gamma0,gammac,Spress,Dzz,sp_energy)
      real(kind=8) :: eta(32),eqSprain,Spress(8,1),etav(32,1),depsdsprain(1,8)
      real(kind=8) :: sprain(8,1), K_sp0(8,8),kappa0,Spress0(8,1),ga,domega2dgamma
      real(kind=8) :: Sprainlocal(2),Dzz(8,8),bmatb(8,32),omega2,gamma0,gammac,sp_energy
      integer:: progDam
      PARAMETER (ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0)
      K_sp0 = 0.0D0
      tmp = 0.0D0
      do i = 1,8
            K_sp0(i,i) = kappa0
      end do 
      etav(:,1) = eta(:) 
      sprain = matmul(bmatb,etav)! sprain 8x1
      call EquivSprain(sprain, eqSprain, depsdsprain)
      Spress0 = matmul(K_sp0,sprain)!*exp(-(eqSprain/0.002d0)) ! spress 8x1
      ga = Sprainlocal(1)
C      gamma = MAX(gamma0,Sprainlocal(1))
c      IF (eqSprain>0.001D0) THEN
c      omega2 =1.0D0
c      else
c      omega2 =1.0D-10
c      END IF
C      write(6,*), "Sprainlocal=",Sprainlocal
C      write(6,*), "eqSprain=",eqSprain
C      write(6,*), "gamma=",ga
      
	IF(eqSprain.GE.ga) THEN
            ga = eqSprain
            progDam =1
      else
            progDam =0
      END IF
      call Damage(ga,omega2,domega2dgamma,gamma0,gammac)
      Sprainlocal(1)=eqSprain!!!!!!!!!!!!!!!!
C      omega2 = MAX(omega2, 1d-12)
c      Sprainlocal(2)=MAX(omega2,Sprainlocal(2))
      Sprainlocal(2)=omega2
      Dzz = Sprainlocal(2)*K_sp0
      IF(progDam == 1) THEN      
            Dzz = Dzz+domega2dgamma*MATMUL(Spress0,depsdsprain)
      END IF
C            write(6,*), "omega2=",omega2
C            write(6,*), "K_sp0=",K_sp0
C            write(6,*), "domega2dgamma=",domega2dgamma
C            write(6,*), "Spress0=",Spress0
C            write(6,*), "depsdsprain=",depsdsprain
C            D=MIN(D,0.99D0)


c      write(6,*), "omega2=",omega2
C      write(6,*), "Dzz=",Dzz
      Spress = omega2*Spress0
      sp_energy = 0.0d0
      DO I =1,8
      sp_energy= sp_energy+0.5D0*Spress(I,1)*sprain(I,1)
      END DO
      end subroutine Spress_Sprain
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine EquivSprain(Sprain, eqSprain, depsdsprain)
      real(kind=8) :: Sprain(8),eqStrain,depsdsprain(1,8),Sprainv(8,1)
      real(kind=8) :: eqSprainm(1,1),eqSprain
      Sprainv(:,1) = Sprain
C      write(6,*), "Sprainv=",Sprainv
      eqSprainm = matmul(transpose(Sprainv),Sprainv)
      eqSprain =eqSprainm(1,1)**0.5d0
      depsdsprain = 0.0D0
      do i = 1,8
         if (eqSprain > 1D-16) then
           depsdsprain(1,i) = Sprain(i)/eqSprain
         end if
      end do

      end subroutine EquivSprain
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      subroutine finalmatrix(K_s, K_f)
      real(kind=8):: K_s(88, 88),K_f(88, 88)
      integer(kind=8):: ind(88)
      ind =   (/1,2,3,25,26,27,28,57,58,59,60,
     $          4,5,6,29,30,31,32,61,62,63,64,
     $          7,8,9,33,34,35,36,65,66,67,68,
     $          10,11,12,37,38,39,40,69,70,71,72,
     $          13,14,15,41,42,43,44,73,74,75,76,
     $          16,17,18,45,46,47,48,77,78,79,80,
     $          19,20,21,49,50,51,52,81,82,83,84,
     $          22,23,24,53,54,55,56,85,86,87,88/)
      do i = 1,88
        do j = 1,88
           K_f(i,j)=K_s(ind(i),ind(j))
        end do
      end do
      
      end subroutine finalmatrix
cccccccccccccccccccccccc
      subroutine Damage(kappa,omega, domegadkappa,kappa0,kappac )
      real(kind=8) :: kappa, kappa0, kappac
      real(kind=8) :: omega, domegadkappa,kappa_d0,kappa_d
      kappa_d0 = 0.4d-3
      kappa_d  = 5d-2
c      write(6,*), "kappa0=",kappa0
c      write(6,*), "kappac=",kappac
c      write(6,*), "kappa=",kappa
c      write(6,*), "omega=",omega
      if (kappa <= kappa0) then
            omega = 1d-12
            domegadkappa = 0.0d0
c            write(6,*), "1"
      elseif (kappa0 < kappa .and. kappa <= (kappac-1.d-13)) then
            fac = kappac / kappa
            omega =1d-12+ fac * (kappa - kappa0) / (kappac -kappa0)
C            write(6,*), "kappa=",kappa
            domegadkappa = fac / (kappac - kappa0) - (omega / kappa)
c            write(6,*), "2"
      else!if ((kappac-1.d-13) < kappa .and. kappa < kappa_d0) then
            omega = 1d-12+(1.0d0-kappa0/((kappac-1.d-13)))/(1-kappa0/kappac)
            domegadkappa = 0.0d0
c            write(6,*), "3"
c      elseif ((kappac-1.d-13) < kappa_d0 .and. kappa < kappa_d-1.d-13) then
c            fac = kappa_d / kappa
c            omega =1.d0 - fac * (kappa - kappa_d0) / (kappa_d -kappa_d0)
c            domegadkappa = -fac / (kappa_d - kappa_d0) + ((1.0d0-omega) / kappa)
c      else 
c            omega = 1d-12
c            domegadkappa = 0.0d0
      end if

      end subroutine Damage
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Damage2(kappa,omega, domegadkappa,kappa0,kappac )
      real(kind=8) :: kappa, kappa0, kappac
      real(kind=8) :: omega, domegadkappa,scale
      REAL, PARAMETER :: PI = 3.14159265358979323846264338327950288419716939937510
c      write(6,*), "kappa0=",kappa0
c      write(6,*), "kappac=",kappac
c      write(6,*), "kappa=",kappa
c      write(6,*), "omega=",omega
      if (kappa <= kappa0) then
            omega = 1d-12
            domegadkappa = 0.0d0
c            write(6,*), "1"
      elseif (kappa0 < kappa .and. kappa < kappac) then
            fac = kappac / kappa
            omega =1d-12+ fac * (kappa - kappa0) / (kappac -kappa0)
C            write(6,*), "kappa=",kappa
            domegadkappa = fac / (kappac - kappa0) - (omega / kappa)

c            write(6,*), "3"
      end if

      end subroutine Damage2