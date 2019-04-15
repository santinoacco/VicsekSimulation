! COMPILAMOS
! \> g95 -Wall -o Vicsek2D ModuloSimulaciones.o ModuloSim_Vicsek.o Vicsek.f90

program VicsekVectorialModel
    use Simulaciones, only: WestRecurence,ACTc_Barkema,SalidaDatosTabla_tvsX
    use Sim_Vicsek
implicit none
integer::ptos,N
integer:: t,tfin,temp,dimMido,midoCada,dimACTc
real:: L,V0,rad,param,eta_0,eta_f,alfa, G
real,allocatable:: r(:,:),v(:,:),tV(:),eta(:),vx(:),vy(:)
character(999):: arch_inParam, arch_outSim, arch_outSim_Sqr, auxread,arch_ACT,auxG
character(2):: choise,corrida
real,allocatable::dv(:,:),phi_mod(:,:),phiSqr_mod(:,:),ACTc(:,:),tt(:)


!======================= ACLARACIONES =======================!
! LxL: dim. de la superficie
! rad: radio de vision de la particula
! V0: modulo de velocidad
! phi: parametro de control
! eta: ruido 'térmico' é [0,1]
! psi: ruido angular é [-Pi,Pi]
!***********************************************************!


! ******* Iniciar la semilla para RANMAR ******* !
!---------------------------------------------------------
integer :: values(1:8), k
integer, dimension(:), allocatable :: seed
call date_and_time(values=values)
call random_seed(size=k)
allocate(seed(1:k))
seed(:)=values(8)
call random_seed(put=seed)
!-------------------- Iniciar RANMAR ---------------------
call RMARIN(seed,0,0)
!************************************************!



! ******* INGRESAR PARAMETROS DEL SISTEMA ******* !
! arch_inParam='Entrada_Parametros_Vicsek.dat'
read(*,*) arch_inParam
open(200,file=arch_inParam)
read(200,*) auxread, N
read(200,*) auxread, L
read(200,*) auxread, V0
read(200,*) auxread, rad
read(200,*) auxread, ptos
read(200,*) auxread, eta_0
read(200,*) auxread, eta_f
read(200,*) auxread, tfin
read(200,*) auxread, midoCada
read(200,*) auxread, corrida
read(200,*) auxread, arch_outSim
! read(*,*) auxread, N
! read(*,*) auxread, V0
! read(*,*) auxread, L
! read(*,*) auxread, rad
! read(*,*) auxread, ptos
! read(*,*) auxread, eta_0
! read(*,*) auxread, eta_f
! read(*,*) auxread, tfin
! read(*,*) auxread, midoCada
! read(*,*) auxread, corrida
! read(*,*) auxread, arch_outSim


close(200)


write(*,*) "Choose a Model between Angular/Vectorial (A/V)"
read(*,*) choise
auxread=trim(choise)//trim(corrida)

select case(choise)
    case("A")
        write(*,*) "Vicsek Angular initialized"
        arch_outSim_Sqr = trim(arch_outSim)//trim(auxread)//"_Sqr.dat"
        arch_outSim = trim(arch_outSim)//trim(auxread)//".dat"
        arch_ACT= "ACorr_V"//trim(auxread)//".dat"
    case("V")
        write(*,*) "provide G"
        read(*,*) G
        write(*,*) "Vicsek Vectorial initialized"

        auxG="_G"//trim(char(floor(G)))//'p'//trim(char(int(G-int(G))*10))
        alfa=ceiling((G-int(G))*10)
        ! write(*,*)alfa
        write(auxG,"(A2,I1,A1,I1)") "_G",int(G),'p',int(alfa)
        write(*,*) auxG
        arch_outSim_Sqr = trim(arch_outSim)//trim(auxread)//trim(auxG)//"_Sqr.dat"
        arch_outSim = trim(arch_outSim)//trim(auxread)//trim(auxG)//".dat"
        arch_ACT= "ACorr_V"//trim(auxread)//trim(auxG)//".dat"
    end select

param=L/2
if (eta_0<0 .or. eta_f>1) then
    stop 'ingreso datos fuera de rango'
end if

! *********************************************** !

allocate(r(N,2))
allocate(v(N,2))

call InitializeState(r,v,V0,L,.false.)


allocate(eta(ptos))
allocate(dv(N,2))

dimMido=int(tfin/midoCada)
allocate(tV(dimMido))
allocate(vx(dimMido))
allocate(vy(dimMido))
! allocate(vvx(dimMido))
! allocate(vvy(dimMido))
allocate(phi_mod(dimMido,ptos))
allocate(phiSqr_mod(dimMido,ptos))



phi_mod = 0
phiSqr_mod = 0
eta=0
write(*,*) 'OutFiles:'
write(*,*) arch_outSim
write(*,*) arch_outSim_Sqr
write(*,*) arch_ACT
write(*,*) '=============================='
write(*,*) 'ruido creciente o decreciente (c/d)?'
read(*,*) auxread

! VARIAMOS EL RUIDO
do temp=1,ptos !cantidad de puntos con distinto ruido
    alfa=(real(temp)-1)/(real(ptos)-1)
    if (auxread == 'c') then
        eta(temp) = eta_0 + alfa*(eta_f-eta_0) !va de eta_0 -> eta_f
        else if (auxread == 'd') then
            eta(temp) = eta_f + alfa*(eta_0-eta_f) !va de eta_f -> eta_0 
        else 
            stop('parametro fuera de rango')
        end if

    write(*,*) 'eta:', eta(temp)
    write(*,*)
    dv=0
    G=1.2
    do t=1,tfin
        
        ! a partir de acá lo que sería un MCS, xq muevo todas las part a la vez
        select case(choise)
            case("A")
                call DinamicaVicsek2D_Angular(r,v,param,rad,L,V0,eta(temp))
            case("V")
                call DinamicaVicsek2D_Vectorial(r,v,param,rad,L,V0,eta(temp),G)
            case default
                stop('please choose a dynamic')
            end select

        if (mod(t,midoCada)==0) then

            j=int(t/midoCada) !indice para guardar la medida
            ! MEDIR cada j cantidad de pasos
            tV(j)=t

            vx=0
            vy=0
            ! vvx=0
            ! vvy=0
            do i=1,N
                vx(j)=vx(j)+v(i,1)
                vy(j)=vy(j)+v(i,2)
                ! vvx(j)=vvx(j) + v(i,1)*v(i,1)
                ! vvy(j)=vvy(j) + v(i,2)*v(i,2)
                enddo

            ! CALCULAR P.O temporal
            phi_mod(j,temp)= sqrt(vx(j)*vx(j)+vy(j)*vy(j))/real(N*V0)
            ! phiSqr_mod(j,temp)= sqrt(vvx(j)*vvx(j)+vvy(j)*vvy(j))/real((N*V0)**2)
            ! phiSqr_mod(j,temp)= sqrt((vx(j)*vx(j)+vy(j)*vy(j))**2)/real((N*V0)**2)
            phiSqr_mod(j,temp)= phi_mod(j,temp)**2 ! esto es lo que creo que corresponde.. consultar por las dudas
            end if
        end do
    end do
deallocate(r)
deallocate(v)
deallocate(vx)
deallocate(vy)
! deallocate(vvx)
! deallocate(vvy)
deallocate(dv)
! ================ CALCULAMOS AutoCorrelacionTemporal ================ !
! dimACTc=dim_vec_eq/2!esto es xq se suele calcular hasta la mitad del vector, luego tengo pocos datos
! dimACTc=dim_vec_eq/2  !checkear q no sea dim_vec_eq/2-1

dimACTc=dimMido/2
allocate(ACTc(dimACTc,ptos))
allocate(tt(dimACTc))
ACTc=0
do temp=1,ptos
    call ACTc_Barkema(phi_mod(:,temp),ACTc(:,temp))
    end do
do t = 1,dimACTc
    tt(t)=t
    enddo

! ESCRIBIR EN EL ARCHIVO DE SALIDA LOS DATOS
open(250,file=arch_outSim)
call SalidaDatosTabla_tvsX(250,tV,phi_mod,eta,'(A3,999(X,f5.2,X))','(A1,999(X,A5,A1,f4.2,A1))','eta','|phi|')
close(250)

open(270,file=arch_outSim_Sqr)
call SalidaDatosTabla_tvsX(270,tV,phiSqr_mod,eta,'(A3,999(X,f5.2,X))','(A1,999(X,A9,A1,f4.2,A1))','eta','|phi^2|')
close(270)

open(300,file=arch_ACT)
call SalidaDatosTabla_tvsX(300,tt,ACTc,eta,'(A3,999(X,f5.2,X))','(A1,999(X,A4,A1,f4.2,A1))','eta','ACTc')
close(300)
deallocate(ACTc)
deallocate(tt)
! =============================== FIN ================================ !
deallocate(tV)
deallocate(eta)
deallocate(phi_mod)
deallocate(phiSqr_mod)



deallocate(seed)


end program


! ****************** SUBRUTINAS AUXILIARES ***************** !
!============================================================!
!============================================================!


