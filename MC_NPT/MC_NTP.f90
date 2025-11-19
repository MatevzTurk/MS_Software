module parametri_simulacije
	integer:: st_paralelk, st_korakov
	real*8:: maks_premik, maks_dV					!maksimalen premik, maksimalna sprememba volumna
	real*8:: Tmin, Tmax, dT							!temperature simulacij
	real*8:: Pmin, Pmax, dP							!tlaki simulacij
end module parametri_simulacije

module box
	integer:: N										!Število delcev
	real*8:: l										!Dolžina škatle
	real*8:: T										!Temperatura v škatli
	real*8:: rho									!Številska gostota delcev
	real*8:: rho0									!začetna gostota delcev
	real*8:: U										!Notranja energija delcev v škatli
	real*8:: P										!tlak v škatli
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
	real*8, allocatable:: rhoavg(:)					!gostota
	real*8, allocatable:: Vavg(:)					!časovna povprečja volumna sistema
	real*8, allocatable:: V2(:)						!časovna povprečja kvadratov volumna sistema
	real*8, allocatable:: Havg(:)					!entalpija
	real*8, allocatable:: H2(:)						!kvadrat entalpije
	real*8, allocatable:: VH(:)						!volumen*entalpija
end module meritve

module simulacija
	integer:: paralelka								!števec paralelk
	integer:: korak									!števec korakov
	character(len=4):: temperatura
	character(len=4):: tlak
	real*8, parameter:: pi=atan(1.0)*4
	integer, parameter:: zacetek_vzorcenja=5000000
	integer, parameter:: korak_vzorcenja=1000
	real*8:: delez_sprejetih_premikov, vsi_premiki
	real*8:: delez_sprejetih_rescaleov, vsi_rescalei
end module simulacija

program MC_NTP

!deklaracija

use parametri_simulacije
use box
use particles
use meritve
use simulacija

implicit none
character:: tekst
integer:: i

open(10, file="MC_NTP_parametri.txt")
read(10, *)
read(10, *)
read(10, *)
read(10, *)
read(10, "(a20, i6)") tekst, N
read(10, "(a20, i6)") tekst, st_paralelk
read(10, "(a20, i10)") tekst, st_korakov
read(10, "(a20, f7.5)") tekst, maks_premik
read(10, "(a20, f7.5)") tekst, maks_dV
read(10, *)
read(10, *)
read(10, *)
read(10, "(a20, f6.3)") tekst, Tmin
read(10, "(a20, f6.3)") tekst, Tmax
read(10, "(a20, f6.3)") tekst, dT
read(10, "(a20, f6.3)") tekst, Pmin
read(10, "(a20, f6.3)") tekst, Pmax
read(10, "(a20, f6.3)") tekst, dP
read(10, "(a20, f6.3)") tekst, rho0
close(10)

T=Tmin
P=Pmin
Tmax=Tmax+dT/2
Pmax=Pmax+dP/2
allocate(delci(N))
do i =1, N
	allocate(delci(i)%pot(N))
end do

allocate(Uavg(st_paralelk))
allocate(U2(st_paralelk))
allocate(rhoavg(st_paralelk))
allocate(Vavg(st_paralelk))
allocate(V2(st_paralelk))
allocate(Havg(st_paralelk))
allocate(H2(st_paralelk))
allocate(VH(st_paralelk))

Uavg=0
U2=0
rhoavg=0
Vavg=0
V2=0
Havg=0
H2=0
VH=0
do i=1, N
	delci(i)%pot=0
end do
delez_sprejetih_premikov=0
vsi_premiki=0
delez_sprejetih_rescaleov=0
vsi_rescalei=0
!simulacija

do while (P<Pmax)
	write(tlak, "(f4.2)") P
	open(100, file="DataP"//tlak//".dat")
	write(100, *) "#Rezultati simulacije Monte Carlo pri tlaku P = "//tlak
	write(100, "(15a15)") "#    T", "rho", "srho", "U", "sU", "H", "sH", "cp", "scp", "k", "sk", "a", "sa", "%x", "%V"
	
	do while (T<Tmax)
		write(temperatura, "(f4.2)") T
		
		do paralelka=1, st_paralelk

			call postavitev()
			
			do korak=1, st_korakov
				call MC_korak()

				if ((korak>zacetek_vzorcenja).and.(mod(korak,korak_vzorcenja)==0)) call sample() 	!vzorčenje
			end do
			
		end do
		call zapis_rezultatov()
		T=T+dT
	end do
	
	close(100)
	P=P+dP
	T=Tmin
	
end do

end program MC_NTP

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
	
	!izračun začetne velikosti škatle
	l=sqrt(real(N)/rho0)
	
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

!premakne naključno izbran delec ali reskalira sistem
subroutine MC_korak()									
	use parametri_simulacije
	use box
	use particles
	use meritve
	use simulacija
	implicit none
	real*8:: a
	integer:: delec									!delec, ki ga premikamo. če out of range, premikamo volumen
	
	call random_number(a)
	delec=int(a*(N+1))+1
	
	if (delec<=N) then								!premik delca
		call premik_delca(delec)
		vsi_premiki=vsi_premiki+1.
		
	else											!premik volumna
		call rescale()
		vsi_rescalei=vsi_rescalei+1.
	end if
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

subroutine rescale()
	use parametri_simulacije
	use box
	use particles
	use meritve
	use simulacija
	implicit none
	real*8:: U_old									
	real*8:: pot_old(N, N)							
	real*8:: l_old, V_old							
	real*8:: dV
	real*8:: V_new
	real*8:: dH
	real*8:: kriterij
	integer:: i
	
	U_old=U
	l_old=l
	V_old=l**2
	do i=1, N
		pot_old(i, :)=delci(i)%pot(:)
	end do
	
	call random_number(dV)
	V_new=V_old+maks_dV*(2.*dV-1.)
	l=sqrt(V_new)

	do i=1, N
		delci(i)%r=delci(i)%r*l/l_old
	end do
	call potenciali()
	
	dH=U-U_old+P*(V_new-V_old)-real(N)*T*log(V_new/V_old)
	call random_number(kriterij)
	if (kriterij>exp(-dH/T)) then
		do i=1, N
			delci(i)%r=delci(i)%r*l_old/l
			delci(i)%pot(:)=pot_old(i, :)
		end do
		U=U_old
		l=l_old
	else
		delez_sprejetih_rescaleov=delez_sprejetih_rescaleov+1.
	end if
end subroutine rescale

!vzorči tekom simulacije
subroutine sample()
	use parametri_simulacije
	use box
	use particles
	use simulacija
	use meritve
	implicit none
	integer:: i, j
	
	Uavg(paralelka)=Uavg(paralelka)+U
	U2(paralelka)=U2(paralelka)+U**2
	rhoavg(paralelka)=rhoavg(paralelka)+real(N)/l**2
	Vavg(paralelka)=Vavg(paralelka)+l**2
	V2(paralelka)=V2(paralelka)+l**4
	Havg(paralelka)=Havg(paralelka)+(U+P*l**2)
	H2(paralelka)=H2(paralelka)+(U+P*l**2)**2
	VH(paralelka)=VH(paralelka)+(U+P*l**2)*l**2
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
	real*8:: H_kon, s_H_kon
	real*8:: rho_kon, s_rho_kon
	real*8:: cp(st_paralelk)
	real*8:: cp_kon, s_cp_kon
	real*8:: k(st_paralelk)
	real*8:: k_kon, s_k_kon
	real*8:: a(st_paralelk)
	real*8:: a_kon, s_a_kon
	real*8:: V_kon, VH_kon, V2_kon, H2_kon
	integer:: i
	
	Uavg=Uavg/real((st_korakov-zacetek_vzorcenja)/korak_vzorcenja)
	U2=U2/real((st_korakov-zacetek_vzorcenja)/korak_vzorcenja)
	rhoavg=rhoavg/real((st_korakov-zacetek_vzorcenja)/korak_vzorcenja)
	Vavg=Vavg/real((st_korakov-zacetek_vzorcenja)/korak_vzorcenja)
	V2=V2/real((st_korakov-zacetek_vzorcenja)/korak_vzorcenja)
	Havg=Havg/real((st_korakov-zacetek_vzorcenja)/korak_vzorcenja)
	H2=H2/real((st_korakov-zacetek_vzorcenja)/korak_vzorcenja)
	VH=VH/real((st_korakov-zacetek_vzorcenja)/korak_vzorcenja)
	
	cp=(H2-Havg**2)/(real(N)*T**2)
	k=(V2-Vavg**2)/(Vavg*T)
	a=(VH-Vavg*Havg)/(Vavg*T**2)
	
	U_kon=0
	U2_kon=0
	H_kon=0
	rho_kon=0
	cp_kon=0
	k_kon=0
	a_kon=0
	do i=1, st_paralelk
		U_kon=U_kon+Uavg(i)
		U2_kon=U2_kon+U2(i)
		H_kon=H_kon+Havg(i)
		rho_kon=rho_kon+rhoavg(i)
		cp_kon=cp_kon+cp(i)
		k_kon=k_kon+k(i)
		a_kon=a_kon+a(i)	
		
		V_kon=V_kon+Vavg(i)
		V2_kon=V2_kon+V2(i)
		VH_kon=VH_kon+VH(i)
		H2_kon=H2_kon+H2(i)
	end do
	U_kon=U_kon/real(st_paralelk)
	U2_kon=U2_kon/real(st_paralelk)
	H_kon=H_kon/real(st_paralelk)
	rho_kon=rho_kon/real(st_paralelk)
	cp_kon=cp_kon/real(st_paralelk)
	k_kon=k_kon/real(st_paralelk)
	a_kon=a_kon/real(st_paralelk)
	
	V_kon=V_kon/real(st_paralelk)
	V2_kon=V2_kon/real(st_paralelk)
	VH_kon=VH_kon/real(st_paralelk)
	H2_kon=H2_kon/real(st_paralelk)
	
	s_U_kon=0
	s_U2_kon=0
	s_H_kon=0
	s_rho_kon=0
	s_cp_kon=0
	s_k_kon=0
	s_a_kon=0
	do i=1, st_paralelk
		s_U_kon=s_U_kon+(Uavg(i)-U_kon)**2
		s_U2_kon=s_U2_kon+(U2(i)-U2_kon)**2
		s_H_kon=s_H_kon+(Havg(i)-H_kon)**2
		s_rho_kon=s_rho_kon+(rhoavg(i)-rho_kon)**2
		s_cp_kon=s_cp_kon+(cp(i)-cp_kon)**2
		s_k_kon=s_k_kon+(k(i)-k_kon)**2
		s_a_kon=s_a_kon+(a(i)-a_kon)**2
	end do
	s_U_kon=sqrt(s_U_kon/real(st_paralelk-1))
	s_U2_kon=sqrt(s_U2_kon/real(st_paralelk-1))
	s_H_kon=sqrt(s_H_kon/real(st_paralelk-1))
	s_rho_kon=sqrt(s_rho_kon/real(st_paralelk-1))
	s_cp_kon=sqrt(s_cp_kon/real(st_paralelk-1))
	s_k_kon=sqrt(s_k_kon/real(st_paralelk-1))
	s_a_kon=sqrt(s_a_kon/real(st_paralelk-1))
	
	!acceptance rate
	delez_sprejetih_premikov=delez_sprejetih_premikov/vsi_premiki
	delez_sprejetih_rescaleov=delez_sprejetih_rescaleov/vsi_rescalei
	
	write(100, "(15f15.6)") T, rho_kon, s_rho_kon, U_kon, s_U_kon, H_kon, s_H_kon, cp_kon, s_cp_kon, k_kon, s_k_kon, a_kon, s_a_kon&
	, delez_sprejetih_premikov, delez_sprejetih_rescaleov
	Uavg=0
	U2=0
	rhoavg=0
	Vavg=0
	V2=0
	Havg=0
	H2=0
	VH=0
	delez_sprejetih_premikov=0
	vsi_premiki=0
	delez_sprejetih_rescaleov=0
	vsi_rescalei=0
	
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
	
	open(300, file="snapshotP"//tlak//"T"//temperatura//".dat")
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
	
	open(300, file="snapshotP"//tlak//"T"//temperatura//".dat")
	do i=1, N
		write(300, "(4f8.4)") delci(i)%r(1), delci(i)%r(2), delci(i)%v(1), delci(i)%v(2)
	end do
	close(300)
end subroutine snapshot2










!justcheck
subroutine Tlakec()
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
	
	write(*,*) vsota/(2.*l**2)+rho*T
end subroutine Tlakec