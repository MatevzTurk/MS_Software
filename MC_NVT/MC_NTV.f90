module parametri_simulacije
	integer:: st_paralelk, st_korakov
	real*8:: maks_premik							!maksimalen premik
	real*8:: Tmin, Tmax, dT							!temperature simulacij
	real*8:: rhomin, rhomax, drho					!gostote simulacij
end module parametri_simulacije

module box
	integer:: N										!Število delcev
	real*8:: l										!Dolžina škatle
	real*8:: dl										!inkrement dolžine škatle pri računanju porazdelitvene funkcije
	real*8:: T										!Temperatura v škatli
	real*8:: rho									!Številska gostota delcev
	real*8:: U										!Notranja energija delcev v škatli
end module box

module particles
	type particle
		real*8:: r(2), v(2), a(2)					!pozicije hitrosti in pospeški v smereh 1:x, 2:y, 3:kot zasuka
		real*8:: zacetna_koordinata(2)				!začetna koordinata delca
		real*8:: koncna_koordinata(2)				!končna koordinata delca brez funkcije minimum image -> za računanje difuzije
		real*8, allocatable:: pot(:)				!potencialne interakcije z delci v škatli
	end type particle
	type(particle), allocatable:: delci(:)			!matrika delcev
end module particles

module meritve
	real*8, allocatable:: Uavg(:)					!časovna povprečja energije 
	real*8, allocatable:: U2(:)						!časovna povprečja kvadratov energije
	real*8, allocatable:: Dif(:)					!difuzijski koeficienti
	real*8, allocatable:: Pavg(:)						!tlaki
	integer, parameter:: porazdelitev_koraki=500
	real*8:: porazdelitev(porazdelitev_koraki)		!porazdelitvena funkcija
end module meritve

module simulacija
	integer:: paralelka								!števec paralelk
	integer:: korak									!števec korakov
	character(len=4):: temperatura
	character(len=4):: gostota
	real*8, parameter:: pi=atan(1.0)*4
	integer, parameter:: zacetek_vzorcenja=5000000
	integer, parameter:: korak_vzorcenja=1000
	real*8:: delez_sprejetih_premikov, vsi_premiki
end module simulacija

program MC_NTV

!deklaracija

use parametri_simulacije
use box
use particles
use meritve
use simulacija

implicit none
character:: tekst
integer:: i

open(10, file="MC_NTV_parametri.txt")
read(10, *)
read(10, *)
read(10, *)
read(10, *)
read(10, "(a20, i6)") tekst, N
read(10, "(a20, i6)") tekst, st_paralelk
read(10, "(a20, i10)") tekst, st_korakov
read(10, "(a20, f7.5)") tekst, maks_premik
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
end do

allocate(Uavg(st_paralelk))
allocate(U2(st_paralelk))
allocate(Dif(st_paralelk))
allocate(Pavg(st_paralelk))

Uavg=0
U2=0
Dif=0
do i=1, N
	delci(i)%pot=0
end do
delez_sprejetih_premikov=0
vsi_premiki=0

!simulacija
do while (T<Tmax)
	write(temperatura, "(f4.2)") T
	open(100, file="DataT"//temperatura//".dat")
	write(100, *) "#Rezultati simulacije Monte Carlo pri temperaturi T = "//temperatura
	write(100, "(11a15)") "#    Rho", "U", "sU", "cv", "scv", "D", "sD", "P", "sP", "t", "%"
	
	do while (rho<rhomax)
		write(gostota, "(f4.2)") rho
		l=sqrt(real(N)/rho)
		dl=(l/2.0)/real(porazdelitev_koraki)
		
		do paralelka=1, st_paralelk
			call postavitev()
			
			do korak=1, st_korakov
				call MC_korak()
				if (korak==zacetek_vzorcenja) then													!začetek vzorčenja
					do i=1, N
						delci(i)%zacetna_koordinata=delci(i)%r
						delci(i)%koncna_koordinata=delci(i)%r
					end do
				end if
				if ((korak>zacetek_vzorcenja).and.(mod(korak,korak_vzorcenja)==0)) call sample() 	!vzorčenje
			end do
			
		end do
		call zapis_rezultatov()
		rho=rho+drho
	end do
	
	close(100)
	T=T+dT
	rho=rhomin
	
end do

end program MC_NTV

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
subroutine interakcIJa(delec1, delec2, pot)
	use particles
	implicit none
	type(particle):: delec1, delec2
	real*8, intent(out):: pot
	real*8:: d
	real*8, parameter:: Ur=2.0
	real*8, parameter:: Rr=1.6
	real*8, parameter:: Ua=1.0
	real*8, parameter:: Ra=2.0
	real*8, parameter:: sa2=0.1
	real*8, parameter:: delta=15.0
	pot=Ur/(1.0+exp(delta*(d(delec1%r, delec2%r)-Rr)))-exp(-((d(delec1%r, delec2%r)-Ra)**2)/(2.0*sa2))+1.0/d(delec1%r, delec2%r)**24
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
	real*8:: z(4)
	real*8:: d
	integer:: i, j
	logical:: prekrivanje
	
	!začetne koordinate
	call random_number(x)
	delci(1)%r=(x-0.5)*l
	i=2
	do while (i<=N)
		 call random_number(x)
		 x=(x-0.5)*l
		 prekrivanje=.true.
		 do j=1, i-1
			if (d(delci(j)%r, x) < 0.7) prekrivanje=.false.
		 end do
		if (prekrivanje) then
			delci(i)%r=x
			i=i+1
		end if
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
	
	U=0
	do i=1, N
		do j=i+1, N
			call interakcIJa(delci(i), delci(j), pot)
			delci(i)%pot(j)=pot
			delci(j)%pot(i)=pot
			U=U+pot
		end do
	end do

end subroutine potenciali

!premakne naključno izbran delec
subroutine MC_korak()									
	use parametri_simulacije
	use box
	use particles
	use meritve
	use simulacija
	implicit none
	real*8:: a
	integer:: delec									!delec, ki ga premikamo

	call random_number(a)							!izbira delca
	delec=int(a*N)+1
	
	call premik_delca(delec)
	vsi_premiki=vsi_premiki+1.

end subroutine MC_korak

subroutine premik_delca(delec)
	use parametri_simulacije
	use box
	use particles
	use meritve
	use simulacija
	implicit none
	integer:: delec
	real*8:: dx(2)									!naključen premik
	real*8:: r_old(2), pot_old(N)					!stari podatki o delcu
	real*8:: r_new(2), dU							!novi podatki o delcu
	real*8:: kriterij
	real*8:: pot
	real*8:: image
	integer:: i
	
	call random_number(dx)
	dx=(dx-0.5)*maks_premik

	r_old=delci(delec)%r
	pot_old(:)=delci(delec)%pot(:)
	
	r_new=delci(delec)%r+dx
	
	delci(delec)%r(1)=image(0.0d0, r_new(1))
	delci(delec)%r(2)=image(0.0d0, r_new(2))
	
	do i=1, N
		if (i/=delec) then
			call interakcIJa(delci(delec), delci(i), pot)
			delci(delec)%pot(i)=pot
		end if 
	end do
	
	dU=sum(delci(delec)%pot(:))-sum(pot_old(:))
	call random_number(kriterij)
	if (kriterij<exp(-dU/T)) then
		delci(delec)%koncna_koordinata=delci(delec)%koncna_koordinata+dx
		do i=1, N
			delci(i)%pot(delec)=delci(delec)%pot(i)
		end do
		U=U+dU
		delez_sprejetih_premikov=delez_sprejetih_premikov+1.
	else
		delci(delec)%r=r_old
		delci(delec)%pot=pot_old
	end if
end subroutine premik_delca

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
	real*8:: d, d2, dpot
	real*8:: TempDif
	integer:: i, j
	
	Uavg(paralelka)=Uavg(paralelka)+U
	U2(paralelka)=U2(paralelka)+U**2
	
	TempDif=0
	do i=1, N
		TempDif=TempDif+d2(delci(i)%koncna_koordinata, delci(i)%zacetna_koordinata)
	end do
	Dif(paralelka)=Dif(paralelka)+TempDif/real(korak/korak_vzorcenja)
	
	do i=1, N
		do j=i+1, N
			Pavg(paralelka)=Pavg(paralelka)-dpot(delci(i), delci(j))*d(delci(i)%r, delci(j)%r)
		end do
	end do

	call porazdelitvena_funkcija()
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
	
	!U, U**2, Dif, P
	Uavg=Uavg/real((st_korakov-zacetek_vzorcenja)/korak_vzorcenja)
	U2=U2/real((st_korakov-zacetek_vzorcenja)/korak_vzorcenja)
	Dif=Dif/real((st_korakov-zacetek_vzorcenja)/korak_vzorcenja*N)
	Pavg=Pavg/(real((st_korakov-zacetek_vzorcenja)/korak_vzorcenja)*2.*l**2)+rho*T
	
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
	porazdelitev=porazdelitev/real(st_paralelk*(st_korakov-zacetek_vzorcenja)/korak_vzorcenja*N)/rho*2 
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
	
	!acceptance rate
	delez_sprejetih_premikov=delez_sprejetih_premikov/vsi_premiki
	
	write(100, "(11f15.6)") rho, U_kon, s_U_kon, cv_kon, s_cv_kon, Dif_kon*1d6, s_Dif_kon*1d6, P_kon, s_P_kon&
	, translacijski_parameter, delez_sprejetih_premikov
	Uavg=0
	U2=0
	Dif=0
	porazdelitev=0
	delez_sprejetih_premikov=0
	vsi_premiki=0
	
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