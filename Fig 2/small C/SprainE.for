      subroutine SprainE(Sprainlocal,PROPS,NPROPS,bmata,bmatb,bmatc,nmatb,NDOFEL,
     1 MLVARX,U,DU,D,N,stress,strain,fd,feta,flambda,K_f,COORDS,sp_energy)
      implicit none

      ! Input arrays and variables
      INTEGER ::      NDOFEL,MLVARX
      real(kind=8) :: bmata(3, 8),bmatb(8, 16),bmatc(4, 8),nmatb(4, 16)
      real(kind=8) :: U(NDOFEL),DU(MLVARX,1),D(3, 3),N(4),stress(3),strain(3)
      real(kind=8) :: PROPS(*),COORDS(2,4)
      integer ::       NPROPS
      ! Output arrays
      real(kind=8):: fd(8,1),feta(16,1),flambda(16,1),K_s(40, 40),K_f(40, 40)

      ! Variables
      real(kind=8) :: dd(8,1), eta(16,1), lambda(16,1),nlambda(4,16)
      real(kind=8) :: dd_new(8,1), eta_new(16,1), lambda_new(16,1),dT_new(1,8)
      real(kind=8) :: dT(1,8),etaT(1,16),nlambdaT(16,4)
      real(kind=8) :: kappa, C_p, A(1,1),A0,l_0

      ! Temporary variables
      real(kind=8) :: temp_matrix1(8, 3),temp,mc_result,h_result
      real(kind=8) :: temp_matrix2(4, 1),temp_matrix3(4,1)
      real(kind=8) :: temp_matrix5(1, 1)
!      real(kind=8) :: temp_matrix6(1, 1)
      real(kind=8) :: temp_matrix7(3, 1)
      real(kind=8) :: temp_vector1(8, 1)

      real(kind=8) :: B2TB2(16,16),B3TB3(8,8),N2TN2(16,16),omega2,Dzz(8,8)
      real(kind=8) :: Sprainlocal(2),eqSprain,gamma0,gammac,Spress(8,1)
      real(kind=8) :: sp_energy
      integer(kind=8):: ind_dof1(8),ind_dof2(16),ind_dof3(16)
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
      ind_dof1 =   (/1,2,11,12,21,22,31,32/)
      ind_dof2 =   (/3,4,5,6,13,14,15,16,23,24,25,26,33,34,35,36/)
      ind_dof3 =   (/7,8,9,10,17,18,19,20,27,28,29,30,37,38,39,40/)
      do i = 1,8
         dd(i,1)=U(ind_dof1(i))
      end do
      dT = transpose(dd)
      do i =1,16
         eta(i,1)=U(ind_dof2(i))
      end do
      etaT = transpose(eta)
      do i =1,16
         lambda(i,1)=U(ind_dof3(i))
      end do

      do i = 1,8
      dd_new(i,1)=U(ind_dof1(i))!+DU(ind_dof1(i),1)
      end do
C      write(6,*), "dd_new=",dd_new
      dT_new = transpose(dd_new)
      do i =1,16
      eta_new(i,1)=U(ind_dof2(i))!+DU(ind_dof2(i),1)
      end do
C      write(6,*), "eta_new=",eta_new
      do i =1,16
      lambda_new(i,1)=U(ind_dof3(i))!+DU(ind_dof3(i),1)
      end do
C      write(6,*), "lambda_new=",lambda_new
      nlambda=0d0
      do i = 1,4
      nlambda(1,1+4*(i-1))=N(i)
      nlambda(2,2+4*(i-1))=N(i)
      nlambda(3,3+4*(i-1))=N(i)
      nlambda(4,4+4*(i-1))=N(i)
      end do
      nlambdaT = transpose(nlambda)


      kappa = 300d0
      l_0   = 5d0

      gamma0 = PROPS(5)
      gammac = PROPS(6)  
C      write(6,*), "gamma0=",gamma0
C      write(6,*), "gammac=",gammac    
      B2TB2 = matmul(transpose(bmatb),bmatb)
      B3TB3 = matmul(transpose(bmatc),bmatc)
      N2TN2 = matmul(transpose(nmatb),nmatb)


C      write(6,*), "A(1,1)=",A(1,1)
C      write(6,*), "A0 new=",A0
      ! Calculate fd
      temp_matrix7(1:3,1) =stress 
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
      ! Block dd (8x8)
!      temp_matrix6 = matmul(nlambda,lambda)
      K_s(1:8, 1:8) = matmul(matmul(transpose(bmata),D),bmata)!+
C      write(6,*), "bmata in sprianE=",bmata
!     $              2d0* temp_matrix6(1,1)*B3TB3
C      write(6,*), "K_s(1:8, 1:8)=", K_s(1:8, 1:8)
C      write(6,*), "0.25 K_s(1:8, 1:8)=", 0.25d0*K_s(1:8, 1:8)
      ! Block deta (8x16)
      K_s(1:8, 9:24) = 0d0
c      write(6,*), "K_s(1:8, 9:24)=", K_s(1:8, 9:24)
      ! Block dlambda (8x4)
      K_s(1:8, 25:40) =  matmul(transpose(bmatc),nlambda)
C      write(6,*), "K_s(1:8, 25:40) dlambda=", 0.25*K_s(1:8, 25:40)
      ! Block etad (16x8)
      K_s(9:24, 1:8) = 0d0
c      write(6,*), "K_s(9:24, 1:8)=", K_s(9:24, 1:8)
      ! Block etaeta (16x16)
      K_s(9:24, 9:24) = l_0*l_0*matmul(transpose(bmatb),matmul(Dzz,bmatb))
C      write(6,*), "0.25 K_s(9:24, 9:24)=", 0.25*K_s(9:24, 9:24)
      ! Block etalambda (16x4)
      K_s(9:24, 25:40) =  -matmul(transpose(nmatb),nlambda)
C      write(6,*), "K_s(9:24, 25:28) etalambda=",0.25*K_s(9:24, 25:40)
      ! Block lambdad (4x8) bmatc
      K_s(25:40, 1:8) = matmul(nlambdaT,bmatc)
C      write(6,*), "K_s(25:40, 1:8)=",0.25*K_s(25:40, 1:8)
      ! Block lambdaeta (4x16)
      K_s(25:40, 9:24) = - matmul(nlambdaT,nmatb)
C      write(6,*), "K_s(25:40, 9:24) lambdaeta=",0.25*K_s(25:40, 9:24)
      ! Block lambdalambda (4x4)
      K_s(25:40, 25:40) = 0.d0
c      write(6,*), "K_s(25:40, 25:40)=",K_s(25:40, 25:40)
C      write(6,*), "K_s=",K_s
      call finalmatrix(K_s,K_f)
c      pause

      end subroutine SprainE

*******************************************************************
      subroutine Spress_Sprain(Sprainlocal,eta,bmatb,eqSprain,kappa0,gamma0,gammac,Spress,Dzz,sp_energy)
      real(kind=8) :: eta(16),eqSprain,Spress(8,1),etav(16,1),depsdsprain(1,8)
      real(kind=8) :: sprain(8,1), K_sp0(8,8),kappa0,Spress0(8,1),ga,domega2dgamma
      real(kind=8) :: Sprainlocal(2),Dzz(8,8),bmatb(8,16),omega2,gamma0,gammac,sp_energy
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
      real(kind=8):: K_s(40, 40),K_f(40, 40)
      integer(kind=8):: ind(40)
      ind =   (/1,2,9,10,11,12,25,26,27,28,
     $          3,4,13,14,15,16,29,30,31,32,
     $          5,6,17,18,19,20,33,34,35,36,
     $          7,8,21,22,23,24,37,38,39,40/)
      do i = 1,40
        do j = 1,40
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