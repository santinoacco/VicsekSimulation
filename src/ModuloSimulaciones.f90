! COMPILAMOS
! \> g95 -Wall -c ModuloSimulaciones.f90

module Simulaciones
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)
  real,parameter:: pi=4.d0*atan(1.d0)
  integer::i,j
  contains
    subroutine SalidaDatosTabla_tvsX(UnitOut,t,Fx,X,HeaderFormat,RowHeaderFormat,XName,FxName)
      integer,intent(in):: UnitOut
      integer::nRow,nCol
      character(*),intent(in):: HeaderFormat,RowHeaderFormat,XName,FxName
      real,intent(in)::t(:),Fx(:,:),X(:)

      nRow = size(t)
      nCol = size(X) !checkear
      write(UnitOut,HeaderFormat) XName,(X(j),j=1,nCol)
      write(UnitOut,RowHeaderFormat) "t",(FxName,"(",X(j),")",j=1,nCol)
      do i=1,nRow
          write(UnitOut,*) t(i),(Fx(i,j),j=1,nCol)
          enddo

      return
      end subroutine SalidaDatosTabla_tvsX
    ! ============================== !
    subroutine WestRecurence(X,mu,S)
      ! ******** ACLARACIONES ******** !
      ! X: cantidad de la cual estimamos <X> y S(X)
      ! mu: estimador de valor medio
      ! S: estimador de varianza
      integer:: cant_muestras
      real,intent(in):: X(:)
      real(dp),intent(out):: mu,S
      real(dp)::Q,R
      mu=0
      S=0
      cant_muestras=size(X)
      do i=1,cant_muestras
        Q=X(i)-mu
        R=Q/i
        mu=mu+R
        S=S+(i-1)*Q*R
        enddo
      S=S/(cant_muestras-1)
      return
      end subroutine WestRecurence
    ! ============================== !
    subroutine Gc_Discreto(red,Xmedio,Gc)
      ! subroutine Gc_Discreto(red,Xmedio,Gc,radio)
      integer, intent(in):: red(:,:), Xmedio
      integer:: r,dimred,param
      real,intent(inout)::Gc(:)
      real:: drr,dx,dy
      integer,allocatable:: N(:)

      ! Gc(r) tiene dim= param, lo mismo q N(r)
      dimred=size(red(:,1))

      param=int(dimred/2)
      allocate(N(param))
      ! allocate(radio(param))
      ! Gc=0
      N=0
      do r=1,param
        ! radio(r)=r
        do i=1,dimred
          do j=i,dimred
            dx=1-i
            dy=1-j
            drr=dx**2+dy**2
            ! ------ FORMA EFICIENTE ------ !
            ! dr=dr+param*floor(drr/param**2)*sign(real(1),dr) !revisar si va con +/-
            dx=dx+param*floor(drr/(param**2))*sign(real(1),dx)
            dy=dy+param*floor(drr/(param**2))*sign(real(1),dy)
            ! dy=dy+param*floor((dy*dy)/param**2)*sign(real(1),dy)
            ! drr=dr**2
            drr=dx**2+dy**2

            ! if (drr<=r*r) then
            if (drr==r*r) then
              N(r)=N(r)+1
              Gc(r)=Gc(r)+red(i,j)*red(1,1)
              end if
            enddo
          enddo
        Gc(r)=real(Gc(r)-Xmedio**2)/real(N(r))
        enddo
      deallocate(N)

      ! do i=1,param
      !   write(*,*)Gc(i)
      !   end do


      return
      end subroutine
    ! ============================== !
    subroutine lectura_matrix(M,unitFile,formato)
      integer::dimN,dimM
      integer,intent(in):: unitFile,M(:,:)
      character(80):: formato
      dimN=size(M(:,1))
      dimM=size(M(1,:))

      if (dimN-int(dimN) /= 0) stop 'ingreso una dimension no entera'
      !escribe una matriz M de dimension=dimN
      if (unitFile == 0) then
          do, i=1,dimN
              write(*,formato) ( M(i,j), j=1,dimM)
              enddo
      else
          do, i=1,dimN
          write(unitFile,formato) ( M(i,j), j=1,dimM)
              enddo
          endif
      return
      end subroutine
    ! ============================== !
    subroutine ACT_c(X,ACTc)
        ! integer:: t, tmax,NA
        integer:: t, tmax
        ! real::aux1,norm,Xmean
        real::Xmean
        real,intent(inout)::X(:)
        real, intent(out):: ACTc(:)
        ! real:: Xaux()

        tmax=size(X)

        ! estimamos la media
        Xmean =0
        do i=1,tmax
            Xmean = Xmean + X(i)
            end do
        Xmean = Xmean/tmax

        ! restamos la media a todo el vector
        do i=1,tmax
            X(i)= X(i)-Xmean
            end do

        ! calculamos ACTc
        do t=1,tmax/2
            do j=t,tmax/2-t
                ! ACTc(j-t+1) = ACTc(j-t+1) + X(t)*(X(j+t)-Xmean)
                ACTc(j-t+1) = ACTc(j-t+1) + X(t)*X(j)
                enddo
            ACTc(t)= ACTc(t)/(tmax/2+1-t)
            end do


        ! NA = size(ACTc) = tmax/2
        ! Normalizo a 1. xq se que el primer elemento es el máximo.
        do i=1,tmax/2
            ACTc(i)=ACTc(i)/ACTc(1)
            end do
        return
        ! do i=1,NA
        !     ACTc(i)=ACTc(i)/(ACTc(1)*(NA-i))
        !     end do
        ! return
        end subroutine ACT_c
    ! ============================== !
end module Simulaciones
!============================================================!
subroutine RANMAR(RVEC,LENV)                                      
! =====================================================================
! Universal random number generator proposed by Marsaglia and Zaman     
! in report FSU-SCRI-87-50                                              
!        modified by F. James, 1988 and 1989, to generate a vector      
!        of pseudorandom numbers RVEC of length LENV, and to put in     
!        the COMMON block everything needed to specify currrent state,  
!        and to add input and output entry points RMARIN, RMARUT.       
!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!!!  Calling sequences for RANMAR:                                  ++ 
!!!      CALL RANMAR (RVEC, LEN)   returns a vector RVEC of LEN     ++ 
!!!                   32-bit random floating point numbers between  ++ 
!!!                   zero and one.                                 ++ 
!!!      CALL RMARIN(I1,N1,N2)   initializes the generator from one ++ 
!!!                   32-bit integer I1, and number counts N1,N2    ++ 
!!!                  (for initializing, set N1=N2=0, but to restart ++ 
!!!                    a previously generated sequence, use values  ++ 
!!!                    output by RMARUT)                            ++ 
!!!      CALL RMARUT(I1,N1,N2)   outputs the value of the original  ++ 
!!!                  seed and the two number counts, to be used     ++ 
!!!                  for restarting by initializing to I1 and       ++ 
!!!                  skipping N1*100000000+N2 numbers.              ++ 
!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      DIMENSION RVEC(*)                                                 
      COMMON/RASET1/U(97),C,I97,J97                                     
      PARAMETER (MODCNS=1000000000)                                     
      SAVE CD, CM, TWOM24, NTOT, NTOT2, IJKL                            
      DATA NTOT,NTOT2,IJKL/-1,0,0/                                      
                                                                       
      IF (NTOT .GE. 0)  GO TO 50                                        
                                                                       
!~         Default initialization. User has called RANMAR without RMARIN. 
      IJKL = 54217137                                                   
      NTOT = 0                                                          
      NTOT2 = 0                                                         
      KALLED = 0                                                        
      GO TO 1                                                           
                                                                       
      ENTRY      RMARIN(IJKLIN, NTOTIN,NTOT2N)                          
!~          Initializing routine for RANMAR, may be called before         
!~          generating pseudorandom numbers with RANMAR. The input        
!~          values should be in the ranges:  0<=IJKLIN<=900 OOO OOO       
!~                                           0<=NTOTIN<=999 999 999       
!~                                           0<=NTOT2N<<999 999 999!      
!~  To get the standard values in Marsaglia's paper, IJKLIN=54217137      
!~                                             NTOTIN,NTOT2N=0            
      IJKL = IJKLIN                                                     
      NTOT = MAX(NTOTIN,0)                                              
      NTOT2= MAX(NTOT2N,0)                                              
      KALLED = 1                                                        
!~           always come here to initialize                               
    1 CONTINUE                                                          
      IJ = IJKL/30082                                                   
      KL = IJKL - 30082*IJ                                              
      I = MOD(IJ/177, 177) + 2                                          
      J = MOD(IJ, 177)     + 2                                          
      K = MOD(KL/169, 178) + 1                                          
      L = MOD(KL, 169)                                                  
      WRITE(6,'(A,I10,2X,2I10)') ' RANMAR INITIALIZED:',IJKL,NTOT,NTOT2 
!~       PRINT '(A,4I10)', '   I,J,K,L= ',I,J,K,L                       
      DO 2 II= 1, 97                                                    
      S = 0.                                                            
      T = .5                                                            
      DO 3 JJ= 1, 24                                                    
         M = MOD(MOD(I*J,179)*K, 179)                                   
         I = J                                                          
         J = K                                                          
         K = M                                                          
         L = MOD(53*L+1, 169)                                           
         IF (MOD(L*M,64) .GE. 32)  S = S+T                              
    3    T = 0.5*T                                                      
    2 U(II) = S                                                         
      TWOM24 = 1.0                                                      
      DO 4 I24= 1, 24                                                   
    4 TWOM24 = 0.5*TWOM24                                               
      C  =   362436.*TWOM24                                             
      CD =  7654321.*TWOM24                                             
      CM = 16777213.*TWOM24                                             
      I97 = 97                                                          
      J97 = 33                                                          
!~ C       Complete initialization by skipping                             
!~ C            (NTOT2*MODCNS + NTOT) random numbers                       
      DO 45 LOOP2= 1, NTOT2+1                                           
      NOW = MODCNS                                                      
      IF (LOOP2 .EQ. NTOT2+1)  NOW=NTOT                                 
      IF (NOW .GT. 0)  THEN                                             
        WRITE(6,'(A,I15)') ' RMARIN SKIPPING OVER ',NOW                 
       DO 40 IDUM = 1, NTOT                                             
       UNI = U(I97)-U(J97)                                              
       IF (UNI .LT. 0.)  UNI=UNI+1.                                     
       U(I97) = UNI                                                     
       I97 = I97-1                                                      
       IF (I97 .EQ. 0)  I97=97                                          
       J97 = J97-1                                                      
       IF (J97 .EQ. 0)  J97=97                                          
       C = C - CD                                                       
       IF (C .LT. 0.)  C=C+CM                                           
   40  CONTINUE                                                         
      ENDIF                                                             
   45 CONTINUE                                                          
      IF (KALLED .EQ. 1)  RETURN                                        
                                                                       
!~           Normal entry to generate LENV random numbers                 
   50 CONTINUE                                                          
      DO 100 IVEC= 1, LENV                                              
      UNI = U(I97)-U(J97)                                               
      IF (UNI .LT. 0.)  UNI=UNI+1.                                      
      U(I97) = UNI                                                      
      I97 = I97-1                                                       
      IF (I97 .EQ. 0)  I97=97                                           
      J97 = J97-1                                                       
      IF (J97 .EQ. 0)  J97=97                                           
      C = C - CD                                                        
      IF (C .LT. 0.)  C=C+CM                                            
      UNI = UNI-C                                                       
      IF (UNI .LT. 0.) UNI=UNI+1.                                       
      RVEC(IVEC) = UNI                                                  
!~              Replace exact zeros by uniform distr. *2**-24             
         IF (UNI .EQ. 0.)  THEN                                         
         ZUNI = TWOM24*U(2)                                             
!~              An exact zero here is very unlikely, but let's be safe.   
         IF (ZUNI .EQ. 0.) ZUNI= TWOM24*TWOM24                          
         RVEC(IVEC) = ZUNI                                              
         ENDIF                                                          
  100 CONTINUE                                                          
      NTOT = NTOT + LENV                                                
         IF (NTOT .GE. MODCNS)  THEN                                    
         NTOT2 = NTOT2 + 1                                              
         NTOT = NTOT - MODCNS                                           
         ENDIF                                                          
      RETURN                                                            
!~            Entry to output current status                              
      ENTRY RMARUT(IJKLUT,NTOTUT,NTOT2T)                                
      IJKLUT = IJKL                                                     
      NTOTUT = NTOT                                                     
      NTOT2T = NTOT2                                                    
      RETURN                                                            
      END                                                        
! ============================== !




    ! ! ============================== !
    ! subroutine AutoCorrelacionTempConect(X,ACTc)
    !   ! ******** ACLARACIONES ******** !
    !   ! ACTc: Auto Correlacion Temporal Conectada
    !   ! X: cantidad a estimar la  ACTc
    !   ! X_mdio: estimador recurente del valor medio de X
    !   ! cant_muestras: cantidad de muestras
    !   ! ****************************** !
    !   ! integer,intent(in)::cant_muestras
    !   integer::i,j,fin,cant_muestras
    !   real,intent(in)::X(:)
    !   real,intent(inout)::ACTc(:)
    !   real::X_mdio,Xs


    !   cant_muestras=size(X)
    !   ! write(*,*) cant_muestras
    !   X_mdio=0
    !   do i=1,cant_muestras
    !     X_mdio=X_mdio+X(i)
    !     enddo
    !   X_mdio=X_mdio/cant_muestras

    !   fin = size(ACTc)
    !   write(*,*) fin
    !   ACTc=0
    !   do i=1,fin
    !     Xs = X(i) - X_mdio
    !     do j=i,fin-i   !está bien hacer la suma en este orden??
    !       ! ACTc(j-i)=ACTc(j-i)+Xs*(X(j)-X_mdio)
    !       ACTc(j-i)=ACTc(j-i)+Xs*(X(j)-X_mdio)/real(fin-i+1)
    !       enddo
    !     enddo

    !    do i=1,fin
    !     ACTc(i)=ACTc(i)/real(fin-i+1)
    !     enddo
    !   return
    !   end subroutine AutoCorrelacionTempConect