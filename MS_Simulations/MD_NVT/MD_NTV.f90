module parametri_simulacije
	integer:: st_paralelk, st_korakov
	real*8:: integracijski_korak					!časovni interval
	real*8:: Tmin, Tmax, dT							!temperature simulacij
	real*8:: rhomin, rhomax, drho					!gostote simulacij
end module parametri_simulacije

module NoseHooverParams
	real*8, dimension(10):: Q		!heat bath masses
	real*8, dimension(10):: vksi	!heat bath velocities
	real*8, dimension(10):: G		!heat bath forces
	real*8, parameter:: Tau=1.0d0	!temperature time constant
	real*8:: Nkft					
end module NoseHooverParams

module box
	integer:: N										!Število delcev
	real*8:: l										!Dolžina škatle
	real*8:: dl										!inkrement dolžine škatle pri računanju porazdelitvene funkcije
	real*8:: T										!Temperatura v škatli
	real*8:: rho									!Številska gostota delcev
	real*8:: P										!Tlak v škatli
	real*8:: U										!Notranja energija delcev v škatli
	real*8:: Wkin									!kinetična energija sistema 1:skupna, 2:translacijska, 3:rotacijska
end module box

module particles
	type particle
		real*8:: r(2), v(2), a(2)					!pozicije hitrosti in pospeški v smereh 1:x, 2:y, 3:kot zasuka
		real*8:: zacetna_koordinata(2)				!začetna koordinata delca
		real*8:: koncna_koordinata(2)				!končna koordinata delca brez funkcije minimum image -> za računanje difuzije
		real*8, allocatable:: pot(:)				!potencialne interakcije z delci v škatli
		real*8, allocatable:: F(:, :)				!sile na delec
	end type particle
	type(particle), allocatable:: delci(:)			!matrika delcev
end module particles

module meritve
	real*8, allocatable:: Uavg(:)					!časovna povprečja energije 
	real*8, allocatable:: U2(:)						!časovna povprečja kvadratov energije
	real*8, allocatable:: Dif(:)					!difuzijski koeficienti
	real*8, allocatable:: Pavg(:)					!tlaki
	integer, parameter:: porazdelitev_koraki=500
	real*8:: porazdelitev(porazdelitev_koraki)		!porazdelitvena funkcija
end module meritve

module simulacija
	integer:: paralelka								!števec paralelk
	integer:: korak									!števec korakov
	character(len=4):: temperatura
	character(len=4):: gostota
	real*8, parameter:: pi=atan(1.0)*4
	integer, parameter:: korak_vzorcenja=100
end module simulacija

program MD_NTV

!deklaracija

use parametri_simulacije
use box
use particles
use meritve
use simulacija
use NoseHooverParams

implicit none
character:: tekst
integer:: i, j

open(10, file="MD_NTV_parametri.txt")
read(10, *)
read(10, *)
read(10, *)
read(10, *)
read(10, "(a20, i6)") tekst, N
read(10, "(a20, i6)") tekst, st_paralelk
read(10, "(a20, i10)") tekst, st_korakov
read(10, "(a20, f7.5)") tekst, integracijski_korak
read(10, *)
read(10, *)
read(10, *)
read(10, "(a20, f6.3)") tekst, Tmin
read(10, "(a20, f6.3)") tekst, Tmax
read(10, "(a20, f6.3)") tekst, dT
read(10, "(a20, f6.3)") tekst, rhomin
read(10, "(a20, f6.3)") tekst, rhomax
read(10, "(a20, f6.3)") tekst, drho
close(10)

T=Tmin
rho=rhomin
Tmax=Tmax+dT/2
rhomax=rhomax+drho/2
allocate(delci(N))
do i =1, N
	allocate(delci(i)%pot(N))
	allocate(delci(i)%F(N, 2))
end do

allocate(Uavg(st_paralelk))
allocate(U2(st_paralelk))
allocate(Dif(st_paralelk))
allocate(Pavg(st_paralelk))

Uavg=0
U2=0
Dif=0
Pavg=0
do i=1, N
	delci(i)%pot=0
	delci(i)%F=0
end do

!simulacija
do while (T<Tmax)
	write(temperatura, "(f4.2)") T
	open(100, file="DataT"//temperatura//".dat")
	write(100, *) "#Rezultati simulacije Molekulske dinamike pri temperaturi T = "//temperatura
	write(100, "(10a15)") "#    Rho", "U", "sU", "cv", "scv", "D", "sD", "P", "sP", "t"
	
	do while (rho<rhomax)
		write(gostota, "(f4.2)") rho
		l=sqrt(real(N)/rho)
		dl=(l/2.0)/real(porazdelitev_koraki)
		
		do paralelka=1, st_paralelk
		
			!NoseHoover parametri
			NkfT=2.0d0*dble(N)*T
			Q(1)=NkfT*Tau**2
			do i=2,10
				Q(i)=T*Tau**2
			enddo
			vksi=0.0d0
			G(1)=0.0d0
			do i=1,9
				G(i+1)=(Q(i)*vksi(i)**2-T)/Q(i+1)
			enddo
			
			call postavitev()
			
			do korak=1, st_korakov									!ekvilibracija
				call MD_korak()
			end do
			
			do i=1, N
				delci(i)%zacetna_koordinata=delci(i)%r
				delci(i)%koncna_koordinata=delci(i)%r
			end do
			
			do korak=1, st_korakov
				call MD_korak()										!začetek vzorčenja
				
				if (mod(korak,korak_vzorcenja)==0) call sample() 	!vzorčenje
			end do
			
		end do
		call zapis_rezultatov()
		rho=rho+drho
	end do
	
	close(100)
	T=T+dT
	rho=rhomin
	
end do

end program MD_NTV

!minimum image koordinate
function image(x1, x2)
	use box
	implicit none
	real*8:: x1, x2
	real*8:: image
	image=(x2-x1)-l*nint((x2-x1)/l)
end function image
!razdalja med delcema i j
function d(r1, r2)
	use box
	implicit none
	real*8:: r1(2), r2(2)
	real*8:: image, d
	d=sqrt(image(r1(1), r2(1))**2+image(r1(2), r2(2))**2)
end function d
!kvadrat razdalje med delcema brez minimum image
function d2(r1, r2)
	implicit none
	real*8:: r1(2), r2(2)
	real*8:: d2
	d2=(r2(1)-r1(1))**2+(r2(2)-r1(2))**2
end function d2

!izračuna interakcije med delcema i j
subroutine interakcIJa(delec1, delec2, pot, F)
	use particles
	implicit none
	type(particle):: delec1, delec2
	real*8, intent(out):: pot, F(2)
	real*8:: dpot
	real*8:: image
	real*8:: d
	real*8, parameter:: Ur=2.0
	real*8, parameter:: Rr=1.6
	real*8, parameter:: Ua=1.0
	real*8, parameter:: Ra=2.0
	real*8, parameter:: sa2=0.1
	real*8, parameter:: delta=15.0
	
	pot=Ur/(1.0+exp(delta*(d(delec1%r, delec2%r)-Rr)))-exp(-((d(delec1%r, delec2%r)-Ra)**2)/(2.0*sa2))+1.0/d(delec1%r, delec2%r)**24
	
	dpot=-Ur*delta*exp(delta*(d(delec1%r, delec2%r)-Rr))/(1.0+exp(delta*(d(delec1%r, delec2%r)-Rr)))**2&
	+(d(delec1%r, delec2%r)-Ra)*exp(-((d(delec1%r, delec2%r)-Ra)**2)/(2.0*sa2))/sa2-24./d(delec1%r, delec2%r)**25

	F(1)=dpot*image(delec1%r(1), delec2%r(1))/d(delec1%r, delec2%r)
	F(2)=dpot*image(delec1%r(2), delec2%r(2))/d(delec1%r, delec2%r)
end subroutine interakcIJa

!izračuna odvod potenciala
function dpot(delec1, delec2)
	use particles
	implicit none
	type(particle):: delec1, delec2
	real*8:: dpot
	real*8:: image
	real*8:: d
	real*8, parameter:: Ur=2.0
	real*8, parameter:: Rr=1.6
	real*8, parameter:: Ua=1.0
	real*8, parameter:: Ra=2.0
	real*8, parameter:: sa2=0.1
	real*8, parameter:: delta=15.0
	dpot=-Ur*delta*exp(delta*(d(delec1%r, delec2%r)-Rr))/(1.0+exp(delta*(d(delec1%r, delec2%r)-Rr)))**2&
	+(d(delec1%r, delec2%r)-Ra)*exp(-((d(delec1%r, delec2%r)-Ra)**2)/(2.0*sa2))/sa2-24./d(delec1%r, delec2%r)**25
end function dpot

!nakljično postavi delce v škatlo
subroutine postavitev()			
	use parametri_simulacije
	use box
	use particles
	use simulacija
	use meritve
	implicit none
	real*8:: x(2)
	real*8:: z(2)
	real*8:: d
	integer:: i, j
	logical:: prekrivanje
	real*8:: v_sum(2)
	
	!začetne koordinate
	call random_number(x)
	delci(1)%r=(x-0.5)*l
	i=2
	do while (i<=N)
		 call random_number(x)
		 x=(x-0.5)*l
		 prekrivanje=.true.
		 do j=1, i-1
			if (d(delci(j)%r, x) < 0.8) prekrivanje=.false.
		 end do
		if (prekrivanje) then
			delci(i)%r=x
			i=i+1
		end if
	end do
	
	!začetne hitrosti
	v_sum=0
	do i=1, N
		call random_number(z)
		delci(i)%v(1)=sqrt(-2*log(z(1)))*cos(2*pi*z(2))*sqrt(T)
		delci(i)%v(2)=sqrt(-2*log(z(1)))*sin(2*pi*z(2))*sqrt(T)
		v_sum=v_sum+delci(i)%v
	end do
	do i=1, N
		delci(i)%v=delci(i)%v-v_sum/real(N)
	end do
	
	call potenciali()
end subroutine postavitev

!izračuna potenciale med delci
subroutine potenciali()
	use parametri_simulacije
	use box
	use particles
	use simulacija
	use meritve
	implicit none
	integer:: i, j
	real*8:: pot
	real*8:: F(2)
	
	U=0
	do i=1, N
		do j=i+1, N
			call interakcIJa(delci(i), delci(j), pot, F)
			delci(i)%pot(j)=pot
			delci(j)%pot(i)=pot
			U=U+pot
			delci(i)%F(j, :)=F(:)
			delci(j)%F(i, :)=-F(:)
		end do
	end do
	
	do i=1, N
		delci(i)%a(1)=sum(delci(i)%F(:, 1))
		delci(i)%a(2)=sum(delci(i)%F(:, 2))
	end do
end subroutine potenciali

subroutine MD_korak()
	use parametri_simulacije
	use box
	use particles
	use meritve
	use simulacija
	implicit none
	real*8:: image
	integer:: i
	
	call Nose_Hoover_thermostat()

	do i=1, N
		delci(i)%v=delci(i)%v+delci(i)%a*integracijski_korak/2.
		
		delci(i)%r=delci(i)%r+delci(i)%v*integracijski_korak
		delci(i)%r(1)=image(0.0d0, delci(i)%r(1))
		delci(i)%r(2)=image(0.0d0, delci(i)%r(2))
		delci(i)%koncna_koordinata=delci(i)%koncna_koordinata+delci(i)%v*integracijski_korak
	end do
	call potenciali()
	
	do i=1, N
		delci(i)%v=delci(i)%v+delci(i)%a*integracijski_korak/2.
	end do
	
	call Nose_Hoover_thermostat()
	!call Berendsen_termostat()
end subroutine MD_korak

!Izračuna kinetično energijo sistema
subroutine Kineticna_energija()
	use parametri_simulacije
	use box
	use particles
	use meritve
	use simulacija
	implicit none
	integer:: i
	
	Wkin=0
	do i=1, N
		Wkin=Wkin+(delci(i)%v(1)**2+delci(i)%v(2)**2)/2.
	end do
end subroutine Kineticna_energija

subroutine Tlak()
	use parametri_simulacije
	use box
	use particles
	use meritve
	use simulacija
	implicit none
	real*8:: d
	real*8:: dpot
	real*8:: vsota
	integer:: i, j
	
	vsota=0
	do i=1, N
		do j=i+1, N
			vsota=vsota-dpot(delci(i), delci(j))*d(delci(i)%r, delci(j)%r)
		end do
	end do
	
	P=vsota/(2.*l**2)+rho*T
end subroutine Tlak

!Berendsenov termostat
subroutine Berendsen_termostat()
	use parametri_simulacije
	use box
	use particles
	use meritve
	use simulacija
	implicit none
	real*8:: lambda_trans
	real*8, parameter:: Tau_trans=0.1
	integer:: i
	
	call kineticna_energija()
	
	lambda_trans=sqrt(1+(integracijski_korak/Tau_trans)*(T/(Wkin/N)-1))
	do i=1, N
		delci(i)%v=delci(i)%v*lambda_trans
	end do
end subroutine Berendsen_termostat

subroutine Nose_Hoover_thermostat()
	use parametri_simulacije
	use box
	use particles
	use meritve
	use simulacija
	use NoseHooverParams
	implicit none
	integer:: i
	real*8:: X, AA
	
	call kineticna_energija()
	
	G(1)=(2.0d0*Wkin-Nkft)/Q(1)
	vksi(10)=vksi(10)+G(10)*integracijski_korak/4.0d0
	do i=1, 9
		X=dexp(-vksi(11-i)*integracijski_korak/8.0d0)
		vksi(10-i)=vksi(10-i)*X
		vksi(10-i)=vksi(10-i)+G(10-i)*integracijski_korak/4.0d0
		vksi(10-i)=vksi(10-i)*X
	end do
	
	AA=dexp(-vksi(1)*integracijski_korak/2.0d0)
	
	do i=1, N
		delci(i)%v=delci(i)%v*AA
	end do
	
	G(1)=(2.0d0*Wkin-Nkft)/Q(1)
	do i=1, 9
		X=dexp(-vksi(i+1)*integracijski_korak/8.0d0)
		vksi(i)=vksi(i)*X
		vksi(i)=vksi(i)+G(i)*integracijski_korak/4.0d0
		vksi(i)=vksi(i)*X
		G(i+1)=(Q(i)*vksi(i)**2-T)/Q(i+1)
	end do
	vksi(10)=vksi(10)+G(10)*integracijski_korak/4.0d0
		
end subroutine Nose_Hoover_thermostat

!prešteje delce na razdaljah po intervalu dr
subroutine porazdelitvena_funkcija()
	use parametri_simulacije
	use box
	use particles
	use simulacija
	use meritve
	implicit none
	real*8:: d
	integer:: i, j, a
	
	do i=1, N-1
		do j=i+1, N
			if (d(delci(i)%r, delci(j)%r)<(l/2+dl)) then
				a=int(d(delci(i)%r, delci(j)%r)/dl)
				porazdelitev(a)=porazdelitev(a)+1.0
			end if
		end do
	end do
end subroutine porazdelitvena_funkcija

!vzorči tekom simulacije
subroutine sample()
	use parametri_simulacije
	use box
	use particles
	use simulacija
	use meritve
	implicit none
	real*8:: d2
	real*8:: TempDif
	integer:: i
	
	Uavg(paralelka)=Uavg(paralelka)+U
	U2(paralelka)=U2(paralelka)+U**2
	
	call Tlak()
	Pavg(paralelka)=Pavg(paralelka)+P
	
	TempDif=0
	do i=1, N
		TempDif=TempDif+d2(delci(i)%koncna_koordinata, delci(i)%zacetna_koordinata)
	end do
	Dif(paralelka)=Dif(paralelka)+TempDif/real(korak/korak_vzorcenja)
	
	call porazdelitvena_funkcija()
	
	call Kineticna_energija()
	write(*,*) Wkin/N
end subroutine sample

!zapiše rezultate poskusa pri temperaturi T in gostoti rho
subroutine zapis_rezultatov()
	use parametri_simulacije
	use box
	use particles
	use simulacija
	use meritve
	implicit none
	real*8:: U_kon, s_U_kon
	real*8:: U2_kon, s_U2_kon
	real*8:: Dif_kon, s_Dif_kon
	real*8:: P_kon, s_P_kon
	real*8:: cv(st_paralelk)
	real*8:: cv_kon, s_cv_kon
	real*8:: translacijski_parameter
	integer:: i
	
	Uavg=Uavg/real(st_korakov/korak_vzorcenja)
	U2=U2/real(st_korakov/korak_vzorcenja)
	Dif=Dif/real(st_korakov/korak_vzorcenja*N)
	Pavg=Pavg/real(st_korakov/korak_vzorcenja)
	
	cv=(U2-Uavg**2)/(real(N)*T**2)
	
	U_kon=0
	U2_kon=0
	Dif_kon=0
	cv_kon=0
	P_kon=0
	do i=1, st_paralelk
		U_kon=U_kon+Uavg(i)
		U2_kon=U2_kon+U2(i)
		Dif_kon=Dif_kon+Dif(i)
		P_kon=P_kon+Pavg(i)
		cv_kon=cv_kon+cv(i)
	end do
	U_kon=U_kon/real(st_paralelk)
	U2_kon=U2_kon/real(st_paralelk)
	Dif_kon=Dif_kon/real(st_paralelk)
	P_kon=P_kon/real(st_paralelk)
	cv_kon=cv_kon/real(st_paralelk)
	
	s_U_kon=0
	s_U2_kon=0
	s_Dif_kon=0
	s_P_kon=0
	s_cv_kon=0
	do i=1, st_paralelk
		s_U_kon=s_U_kon+(Uavg(i)-U_kon)**2
		s_U2_kon=s_U2_kon+(U2(i)-U2_kon)**2
		s_Dif_kon=s_Dif_kon+(Dif(i)-Dif_kon)**2
		s_P_kon=s_P_kon+(Pavg(i)-P_kon)**2
		s_cv_kon=s_cv_kon+(cv(i)-cv_kon)**2
	end do
	s_U_kon=sqrt(s_U_kon/real(st_paralelk-1))
	s_U2_kon=sqrt(s_U2_kon/real(st_paralelk-1))
	s_Dif_kon=sqrt(s_Dif_kon/real(st_paralelk-1))
	s_P_kon=sqrt(s_P_kon/real(st_paralelk-1))
	s_cv_kon=sqrt(s_cv_kon/real(st_paralelk-1))
	
	!porazdelitvena funkcija
	porazdelitev=porazdelitev/real(st_paralelk*st_korakov/korak_vzorcenja*N)/rho*2 
	do i=1, porazdelitev_koraki
		porazdelitev(i)=porazdelitev(i)/(pi*((i*dl)**2-((i-1)*dl)**2))
	end do
	
	open(200, file="porazdelitevT"//temperatura//"G"//gostota//".dat")
	do i=1, porazdelitev_koraki
		write(200, *) i*dl, porazdelitev(i)
	end do
	close(200)
	
	translacijski_parameter=0
	do i=1, porazdelitev_koraki-1
		translacijski_parameter=translacijski_parameter+abs(porazdelitev(i)+porazdelitev(i+1)-2.)&
		*sqrt(real(N))/(2.*real(porazdelitev_koraki))/2.
	end do
	
	write(100, "(10f15.6)") rho, U_kon, s_U_kon, cv_kon, s_cv_kon, Dif_kon*1d6, s_Dif_kon*1d6, P_kon, s_P_kon, translacijski_parameter
	Uavg=0
	U2=0
	Dif=0
	Pavg=0
	porazdelitev=0
	
	call snapshot()
end subroutine zapis_rezultatov

!shrani snapshot sistema (x, y, kot)
subroutine snapshot()
	use parametri_simulacije
	use box
	use particles
	use simulacija
	use meritve
	implicit none
	integer:: i
	
	open(300, file="snapshotT"//temperatura//"G"//gostota//".dat")
	do i=1, N
		write(300, "(3f8.4)") delci(i)%r
	end do
	close(300)
end subroutine snapshot

!shrani snapshot sistema (x, y, vx, vy)
subroutine snapshot2()
	use parametri_simulacije
	use box
	use particles
	use simulacija
	use meritve
	implicit none
	integer:: i
	
	open(300, file="snapshotT"//temperatura//"G"//gostota//".dat")
	do i=1, N
		write(300, "(4f8.4)") delci(i)%r(1), delci(i)%r(2), delci(i)%v(1), delci(i)%v(2)
	end do
	close(300)
end subroutine snapshot2