! COMPILAMOS
! \> g95 -Wall -o AnalisisACTc ModuloSimulaciones.o analisisACTc.f90
program pruebaACT
    use Simulaciones, only: ACT_c,SalidaDatosTabla_tvsX
implicit none
integer::ptos,N,i,j
integer:: t,tfin,temp,dimMido,midoCada,dimEq,t_eq
real:: eta_0,eta_f,alfa,tt
real,allocatable:: tV(:),eta(:)
character(999)::  arch_outSim, auxread,arch_ACT,arch_inParam
! character(2)::corrida
real,allocatable::phi_mod(:,:),ACTc(:,:),phi_mod_Eq(:,:)


write(*,*) 'ingrese archivo de parametros'
read(*,*) arch_inParam
open(200,file=arch_inParam)
read(200,*) !auxread, N
read(200,*) !auxread, L
read(200,*) !auxread, V0
read(200,*) !auxread, rad
read(200,*) auxread, ptos
read(200,*) auxread, eta_0
read(200,*) auxread, eta_f
read(200,*) auxread, tfin
read(200,*) auxread, midoCada
! read(200,*) auxread, corrida
read(200,*) auxread, arch_outSim

! ptos=10
! midoCada=200
dimMido = tfin/midoCada
! eta_0=0.1
! eta_f=0.8
write(*,*) 'ingresar archivo a leer'
read(*,*) arch_outSim
write(*,*) 'ingresar archivo salida'
read(*,*) arch_ACT

write(*,*) 'X must be a stationary array'

allocate(eta(ptos))
do temp=1,ptos !cantidad de puntos con distinto ruido
    alfa=(real(temp)-1)/(real(ptos)-1)
    eta(temp) = eta_f + alfa*(eta_0-eta_f) !va de eta_f -> eta_0 
enddo


allocate(tV(dimMido))
allocate(phi_mod(dimMido,ptos))

open(200,file=arch_outSim) 

read(200,*)
read(200,*)
! leer phi

do i=1,dimMido
    read(200,*) tV(i),(phi_mod(i,j),j=1,ptos)
    end do
close(200)

!asumir un tiempo de equilibracion para phi_mod
write(*,*) 'ingrese tiempo de equilibracion'
read(*,*) tt
t_eq=int(tt/midoCada)


dimEq=dimMido-t_eq

! calcular ACTc
allocate(phi_mod_Eq(dimEq,ptos))
do j=1,ptos
    do i=1,dimEq
        phi_mod_Eq(i,j)=phi_mod(i+t_eq,j)
        enddo
    enddo
deallocate(phi_mod)
allocate(ACTc(dimEq,ptos))
! allocate(ACTc(dimMido,ptos))

ACTc=0
do temp=1,ptos
    call ACT_c(phi_mod_Eq(:,temp),ACTc(:,temp))
    end do

deallocate(phi_mod_Eq)


open(300,file=arch_ACT)
! call SalidaDatosTabla_tvsX(300,tV(t_eq:),ACTc,eta,'(A3,999(X,f5.2,X))','(A1,999(X,A4,A1,f4.2,A1))','eta','ACTc')
call SalidaDatosTabla_tvsX(300,tV(:dimEq),ACTc,eta,'(A3,999(X,f5.2,X))','(A1,999(X,A4,A1,f4.2,A1))','eta','ACTc')

close(300)


deallocate(eta)
deallocate(ACTc)
deallocate(tV)


end program
