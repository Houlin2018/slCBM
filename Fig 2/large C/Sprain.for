      subroutine Sprain(bmata,bmatb,bmatc,nmatb,NDOFEL,MLVARX,U,DU,D,N,stress,fd,feta,flambda,K_f)
      implicit none

      ! Input arrays and variables
      INTEGER ::      NDOFEL,MLVARX
      real(kind=8) :: bmata(3, 8),bmatb(8, 16),bmatc(4, 8),nmatb(4, 16)
      real(kind=8) :: U(NDOFEL),DU(MLVARX,1),D(3, 3),N(4),stress(3)

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
      real(kind=8) :: temp_matrix4(16, 16),temp_matrix5(1, 1)
!      real(kind=8) :: temp_matrix6(1, 1)
      real(kind=8) :: temp_matrix7(3, 1)
      real(kind=8) :: temp_vector1(8, 1)

      real(kind=8) :: B2TB2(16,16),B3TB3(8,8),N2TN2(16,16)
      integer(kind=8):: ind_dof1(8),ind_dof2(16),ind_dof3(16)
      integer(kind=8):: i
c      write(6,*), "N shape function=",N	
c      write(6,*), "bmata=",bmata	
c      write(6,*), "bmatb=",bmatb
c      write(6,*), "bmatc=",bmatc
c      write(6,*), "nmatb=",nmatb
c      write(6,*), "U=",U
c      write(6,*), "DU=",DU(:,1)
!      write(6,*), "DU IS",KIND(DU)
!      PAUSE
      fd = 0d0
      feta=0d0
      flambda=0d0
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
      dd_new(i,1)=U(ind_dof1(i))+DU(ind_dof1(i),1)
      end do

      dT_new = transpose(dd_new)
      do i =1,16
      eta_new(i,1)=U(ind_dof2(i))+DU(ind_dof2(i),1)
      end do
      do i =1,16
      lambda_new(i,1)=U(ind_dof3(i))+DU(ind_dof3(i),1)
      end do
      nlambda=0d0
      do i = 1,4
      nlambda(1,1+4*(i-1))=N(i)
      nlambda(2,2+4*(i-1))=N(i)
      nlambda(3,3+4*(i-1))=N(i)
      nlambda(4,4+4*(i-1))=N(i)
      end do
      nlambdaT = transpose(nlambda)

      kappa = 150d0
      l_0   = 25d0
      C_p   = 0.d0

      B2TB2 = matmul(transpose(bmatb),bmatb)
      B3TB3 = matmul(transpose(bmatc),bmatc)
      N2TN2 = matmul(transpose(nmatb),nmatb)

      A = matmul(matmul(transpose(eta_new),B2TB2),eta_new)
c      write(6,*), "eta_new=",eta_new
c      write(6,*), "lambda_new=",lambda_new
c      write(6,*), "B2TB2=",B2TB2
      A0 = A(1,1)**0.5d0
c      write(6,*), "A(1,1)=",A(1,1)
c      write(6,*), "A0 new=",A0
      ! Calculate fd
      temp_matrix7(1:3,1) =stress 
      fd = matmul(transpose(bmata), temp_matrix7)
      fd = fd + matmul(matmul(transpose(bmatc), nlambda),lambda_new)
c      write(6,*), "fd=",fd	

      ! Calculate feta
      call macaulay_brackets(l_0*A0-C_p,mc_result)
c      write(6,*), "l_0*A0-C_p=",l_0*A0-C_p	
      if (mc_result>0) then
      feta = matmul( B2TB2, eta_new)*kappa*mc_result/A0
      else
      feta = 0.0d0
      end if
      feta =feta - matmul(transpose(nmatb),matmul(nlambda,lambda_new))
c      write(6,*), "feta=",feta	
      ! Calculate flambda
      temp_matrix2 = matmul(bmatc,dd_new)
      flambda = matmul(nlambdaT, temp_matrix2)
      temp_matrix3 = matmul(nmatb,eta_new)
      flambda = flambda - matmul( nlambdaT,temp_matrix3)
c      flambda = 0d0
c      write(6,*), "flambda=",flambda
      ! Calculate K_s matrix blocks
      ! K_s matrix blocks:
      ! dd, deta, dlambda
      ! etad, etaeta, etalambda
      ! lambdad, lambdaeta, lambdalambda
      K_s = 0d0
      ! Block dd (8x8)
!      temp_matrix6 = matmul(nlambda,lambda)
      K_s(1:8, 1:8) = matmul(matmul(transpose(bmata),D),bmata)!+
!     $              2d0* temp_matrix6(1,1)*B3TB3
c      write(6,*), "K_s(1:8, 1:8)=", K_s(1:8, 1:8)
      ! Block deta (8x16)
      K_s(1:8, 9:24) = 0d0
c      write(6,*), "K_s(1:8, 9:24)=", K_s(1:8, 9:24)
      ! Block dlambda (8x4)
      K_s(1:8, 25:40) = matmul(transpose(bmatc),nlambda)
c      write(6,*), "K_s(1:8, 25:40)=", K_s(1:8, 25:40)
      ! Block etad (16x8)
      K_s(9:24, 1:8) = 0d0
c      write(6,*), "K_s(9:24, 1:8)=", K_s(9:24, 1:8)
      ! Block etaeta (16x16)
      A = matmul(matmul(transpose(eta),B2TB2),eta)
      A0 = A(1,1)**0.5d0
c      write(6,*), "B2TB2=",B2TB2
c      write(6,*), "A0 old=",A0
c      write(6,*), "eta old=",eta
      call Heaviside(l_0*A0-C_p, h_result)
      temp_matrix4 = matmul(matmul(matmul(B2TB2,eta),etaT),B2TB2)
c      temp_matrix5 = matmul(nlambda,lambda)
      if (h_result > 0) then
      K_s(9:24, 9:24) = kappa*A0**(-2d0)*(l_0*h_result-mc_result/A0)*temp_matrix4+
     $               kappa*A0**(-1d0)*mc_result*B2TB2
      else
      K_s(9:24, 9:24) = 0d0
      end if
c      write(6,*), "K_s(9:24, 9:24)=", K_s(9:24, 9:24)
      ! Block etalambda (16x4)
      K_s(9:24, 25:40) = - matmul(transpose(nmatb),nlambda)
c      write(6,*), "K_s(9:24, 25:28)=", K_s(9:24, 25:40)
      ! Block lambdad (4x8) bmatc
      K_s(25:40, 1:8) = matmul(nlambdaT,bmatc)
c      write(6,*), "K_s(25:40, 1:8)=",K_s(25:40, 1:8)
      ! Block lambdaeta (4x16)
      K_s(25:40, 9:24) = - matmul(nlambdaT,nmatb)
c      write(6,*), "K_s(25:40, 9:24)=",K_s(25:40, 9:24)
      ! Block lambdalambda (4x4)
      K_s(25:40, 25:40) = 0.d0
c     write(6,*), "K_s(25:40, 25:40)=",K_s(25:40, 25:40)
      call finalmatrix(K_s,K_f)
c      pause

      end subroutine Sprain

*******************************************************************
      subroutine macaulay_brackets(x, result)
        implicit none
        real(8), intent(in) :: x
        real(8), intent(out) :: result

        ! Define the Macaulay brackets function
        if (x <= 0.0d0) then
          result = 0.0d0
        else
          result = x
        end if
      end subroutine macaulay_brackets
*******************************************************************
      subroutine Heaviside(x, result)
        implicit none
        real(8), intent(in) :: x
        real(8), intent(out) :: result

        ! Define the Heaviside function
        if (x <= 0.0d0) then
          result = 0.0d0
        else
          result = 1.0d0
        end if
      end subroutine Heaviside
********************************************************************
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