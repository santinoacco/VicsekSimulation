! COMPILAMOS
! \> g95 -Wall -c ModuloSim_Vicsek.f90

module Sim_Vicsek
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer::i,j
  contains
    subroutine InitializeState(r,v,V0,L,printState)
      ! Inicia los vectores 'r' y 'v' aleatoriamente a partir de
      ! los parametros V0 y L. printState es una variable para
      ! decidir si imprimir en pantalla dichos vectores creados
      integer:: N
      real,intent(in):: V0,L
      logical,intent(in):: printState
      real,intent(inout):: r(:,:),v(:,:)
      real, allocatable:: r1(:),r2(:)

      N=size(r(:,1))

      ! ******* GENERAR ESTADO INICIAL ******* !
      allocate(r1(N))
      allocate(r2(N))
      ! --------POSICION-------- !
      ! allocate(r(N,2))
      call ranmar(r1,N)
      call ranmar(r2,N)

      do i=1,N
          r(i,1)=r1(i)*L
          r(i,2)=r2(i)*L
      enddo
      ! --------VELOCIDAD-------- !
      ! allocate(v(N,2))
      do i=1,N
          v(i,1)= V0*cos(r1(i))
          v(i,2)= V0*sin(r1(i))
      enddo
      deallocate(r1)
      deallocate(r2)


      if (printState) then
          write(*,*) 'vector posicion inicial'
          do i=1,N
          write(*,*) r(i,:)
              enddo

          write(*,*) 'vector velocidad inicial'
          do i=1,N
          write(*,*) v(i,:)
              enddo
          end if

      write(*,*) '**************************'
      ! ************************************** !
      return
      end subroutine InitializeState
    ! ============================== !
    elemental real function F1(V0,Vx,Vy,G) result(res)
      real,intent(in)::V0,G,Vx,Vy
      ! real,intent(out)::res
      !G es una constante
      res=1+G*(V0/sqrt(Vx**2+Vy**2)-1)
      end function
    ! ============================== !
    elemental subroutine Distancia(param,X1,X2,dXX)
      !ingreso escalares, devuelve elemento de vector escalar
      real,intent(in):: param,X1,X2
      real,intent(inout):: dXX
      real::dX
      dX=X1-X2 !calcula diferencia
      dXX=dX**2 !calcula dist^2
      ! ------ FORMA EFICIENTE ------ !
      dX=dX+param*floor(dXX/param**2)*sign(real(1),dX) !revisar si va con +/-
      dXX=dX**2
      end subroutine Distancia
    ! ============================== !
    subroutine SumNeighboursVelocity(r,v,dv,param,rad,NNbrs)
      integer::N
      real::dxx,dyy
      real,intent(inout)::r(:,:),v(:,:),dv(:,:),NNbrs(:)
      real,intent(in)::param,rad


      N=size(r(:,1))
      NNbrs=0
      dv=0
      do i=1,N
        do j=i,N
        ! do j=i+1,N !si uso esto tengo que guardar la info tambien en |j-i|
        ! el 'do' anterior no cuenta la particula estudiada, si esta no vió vecinos dv y NNbrs serán =0
        ! esto puede traer errores al dividir por |dv| o por NNbrs.
          ! CALCULO LA DISTANCIA ENTRE TODOS LOS VECINOS
          call Distancia(param,r(i,1),r(j,1),dxx)
          call Distancia(param,r(i,2),r(j,2),dyy)
          ! * notar que la subrutina "distancia" me da dxx como (xi-xj)^2 *

          ! BUSCO TODOS LOS VECINOS EN UN RADIO 'rad' Y SUMO VELOCIDADES
          if (dxx+dyy<rad**2) then
            dv(i,1)= dv(i,1) + v(j,1)  
            dv(i,2)= dv(i,2) + v(j,2)
            dv(j,1)= dv(j,1) + v(i,1)
            dv(j,2)= dv(j,2) + v(i,2)

            NNbrs(i)= NNbrs(i) + 1
            NNbrs(j)= NNbrs(j) + 1
            end if
          end do
        end do
      !DEVUELVE LAS COMPONENTES DE UN VECTOR SIN NORMALIZAR
      return
      end subroutine
    ! ============================== !
    subroutine DinamicaVicsek2D_Angular(r,v,param,rad,L,V0,eta)
      use Simulaciones
      integer::N
      real::V0
      real,allocatable::psi(:),NNbrs(:),dtheta(:),dv(:,:)
      real,intent(inout)::r(:,:),v(:,:)
      real,intent(in)::param,rad,L,eta


      N=size(r(:,1))

      allocate(psi(N))
      allocate(dtheta(N))
      allocate(dv(N,2))
      allocate(NNbrs(N))
      call ranmar(psi,N) !como es una subroutina que llama a cada paso, para cada tiempo se genera una nueva secuencia de ruido
      psi = pi*(2*psi-1) !transforma al vector aleatorio para que esté é [-Pi,Pi]

      dv=0
      
      ! CONTAMOS LOS VECINOS Y SUMAMOS VELOCIDADES
      call SumNeighboursVelocity(r,v,dv,param,rad,NNbrs)
      
      dtheta=0 !solo acá reinicio este contador

      ! CALCULAMOS EL ANGULO DE DESVIACION Y LE AGREGAMOS RUIDO
      do i=1,N
        dtheta(i) = atan2(v(i,2)-dv(i,2),v(i,1)-dv(i,1)) + eta*psi(i)

        ! ACTUALIZO VELOCIDADES
        v(i,1)=V0*cos(dtheta(i))   
        v(i,2)=V0*sin(dtheta(i))  
        end do


      ! ACTUALIZO POSICIÓN DE TODAS LAS PARTICULAS, eligiendo dt=1
      do i=1,N
        r(i,1)= r(i,1)+v(i,1)
        r(i,2)= r(i,2)+v(i,2)
 
        ! CONTROL que no se vayan fuera de la caja

        ! ------ FORMA EFICIENTE ------ !
        r(i,1)=r(i,1)-L*floor(r(i,1)/L)
        r(i,2)=r(i,2)-L*floor(r(i,2)/L)

        end do

      deallocate(psi) 
      deallocate(NNbrs)
      deallocate(dtheta)
      deallocate(dv)

      return
      end subroutine DinamicaVicsek2D_Angular
    ! ============================== !
    subroutine DinamicaVicsek2D_Vectorial(r,v,param,rad,L,V0,eta,G)
      use Simulaciones
      integer::N
      real::V0,G
      real,allocatable::psix(:),psiy(:),denom(:),NNbrs(:),dv(:,:)
      real,intent(inout)::r(:,:),v(:,:)
      real,intent(in)::param,rad,L,eta


      N=size(r(:,1))

      allocate(NNbrs(N))
      allocate(denom(N))
      allocate(psix(N))
      allocate(psiy(N))
      allocate(dv(N,2))

      call ranmar(psix,N)

      psix = pi*(2*psix-1) !transforma al vector aleatorio para que esté é [-Pi,Pi]
      dv=0
      NNbrs=0

      call SumNeighboursVelocity(r,v,dv,param,rad,NNbrs)
      ! COTA DE dv =< V0*NNbrs
      ! write(*,*) 'v antes',v(1,:)
      do i=1,N

        ! LLEVAMOS dv -> dv =< NNbrs
        dv(i,1)=dv(i,1)/real(V0)
        dv(i,2)=dv(i,2)/real(V0)

        ! INTRODUCIMOS ERROR DELTA CORRELACIONADO É [0:1]
        psiy(i) = sin(psix(i))
        psix(i) = cos(psix(i))
        
        !CORREGIMOS EL dv SEGUN Ginelli
        dv(i,1)=dv(i,1) + psix(i)*NNbrs(i)*eta
        dv(i,2)=dv(i,2) + psiy(i)*NNbrs(i)*eta

        ! /real(NNbrs(i))

        denom(i) = sqrt(dv(i,1)**2+dv(i,2)**2)


        ! ACTUALIZO VELOCIDADES, NORMALIZANDO Y LLEVANDO dv -> dv=<V0

        ! ! ---------- esto anda ----------
        ! v(i,1)= V0*dv(i,1)/real(denom(i))
        ! v(i,2)= V0*dv(i,2)/real(denom(i))
        ! ! -------------------------------
        
        v(i,1)= F1(V0,dv(i,1),dv(i,2),G)*dv(i,1)/real(denom(i))
        v(i,2)= F1(V0,dv(i,1),dv(i,2),G)*dv(i,2)/real(denom(i))
        

        end do
        ! write(*,*) 'v despues',v(1,:)
      do i=1,N ! ACTUALIZO LA POSICIÓN DE TODAS LAS PARTICULAS, eligiendo dt=1
        r(i,1)= r(i,1)+v(i,1)
        r(i,2)= r(i,2)+v(i,2)
 
        ! controlas que no se vayan fuera de la caja

        ! ------ FORMA EFICIENTE ------ !
        r(i,1)=r(i,1)-L*floor(r(i,1)/L)
        r(i,2)=r(i,2)-L*floor(r(i,2)/L)

        end do
      deallocate(NNbrs)
      deallocate(denom)
      deallocate(psix)
      deallocate(psiy)
      deallocate(dv)
      return
      end subroutine DinamicaVicsek2D_Vectorial
    ! ============================== !
end module Sim_Vicsek