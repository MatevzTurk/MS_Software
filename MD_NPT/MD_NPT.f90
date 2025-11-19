module parametri_simulacije
	integer:: st_paralelk, st_korakov
	real*8:: integracijski_korak					!časovni interval
	real*8:: Tmin, Tmax, dT							!temperature simulacij
	real*8:: Pmin, Pmax, dP							!tlaki simulacij
end module parametri_simulacije

module box
	integer:: N										!Število delcev
	real*8:: l										!Dolžina škatle
	real*8:: dl										!inkrement dolžine škatle pri računanju porazdelitvene funkcije
	real*8:: T										!Temperatura v škatli
	real*8:: rho									!Številska gostota delcev
	real*8:: rho0									!začetna številska gostota delcev
	real*8:: P										!Nastavljen tlak
	real*8:: P_real									!Pravi tlak v škatli
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
	integer, parameter:: korak_vzorcenja=100
end module simulacija

module NoseHooverParams
	real*8, dimension(10):: Q		!heat bath masses
	real*8, dimension(10):: vksi	!heat bath velocities
	real*8, dimension(10):: G		!heat bath forces
	real*8, parameter:: Tau=0.2d0	!temperature time constant
	real*8:: Nkft					
end module NoseHooverParams

module MTTKParams
	real*8:: Nu, Eps, l0
	real*8:: Qbar
	real*8, parameter:: TauBar=4.0d0
end module MTTKParams

program MD_NTP

!deklaracija

use parametri_simulacije
use box
use particles
use meritve
use simulacija
use NoseHooverParams
use MTTKParams

implicit none
character:: tekst
integer:: i

open(10, file="MD_NPT_parametri.txt")
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
read(10, "(a20, f6.3)") tekst, Pmin
read(10, "(a20, f6.3)") tekst, Pmax
read(10, "(a20, f6.3)") tekst, dP
read(10, "(a20, f6.3)") tekst, rho0
close(10)

T=Tmax
P=Pmin
Tmin=Tmin-dT/2
Pmax=Pmax+dP/2
allocate(delci(N))
do i =1, N
	allocate(delci(i)%pot(N))
	allocate(delci(i)%F(N, 2))
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
	delci(i)%F=0
end do



!simulacija
do while (P<Pmax)
	write(tlak, "(f4.2)") P
	open(100, file="DataP"//tlak//".dat")
	write(100, *) "#Rezultati simulacije Molekulske dinamike pri tlaku P = "//tlak
	write(100, "(13a15)") "#    T", "rho", "srho", "U", "sU", "H", "sH", "cp", "scp", "k", "sk", "a", "sa"
	
	call postavitev()
	Nu=0.0d0
	Eps=0.0d0
	
	
	do while (T>Tmin)
		write(temperatura, "(f4.2)") T
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
			!MTTK parametri
			Qbar=NkfT*TauBar**2
				
		write(*,*) T, P, "init"
		
		do korak=1, st_korakov
			call MD_korak()
			if (mod(korak,korak_vzorcenja)==0) then
				call tlakec()
				call kineticna_energija()
				write(*,*) Wkin/N, P_real
			end if
		end do
		
		do paralelka=1, st_paralelk
		
		write(*,*) T, P, paralelka
		
			do korak=1, st_korakov
				call MD_korak()				
				if (mod(korak,korak_vzorcenja)==0) call sample() 	!vzorčenje
			end do
			
		end do
		call zapis_rezultatov()
		T=T-dT
	end do
	
	close(100)
	P=P+dP
	T=Tmax
end do

end program MD_NTP

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
	use MTTKParams
	implicit none
	real*8:: x(2)
	real*8:: z(2)
	real*8:: d
	integer:: i, j
	logical:: prekrivanje
	real*8:: v_sum(2)
		
	!izračun začetne velikosti škatle
	l=sqrt(real(N)/rho0)
	l0=l
	
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
	use MTTKParams
	implicit none
	real*8:: image
	real*8:: scalar
	integer:: i
	
	call Nose_Hoover_thermostat()
	do i=1, N
		delci(i)%v=delci(i)%v+delci(i)%a*integracijski_korak/2.0d0
	end do
	
	do i=1, N
		delci(i)%r=delci(i)%r+delci(i)%v*integracijski_korak
		delci(i)%r(1)=image(0.0d0, delci(i)%r(1))
		delci(i)%r(2)=image(0.0d0, delci(i)%r(2))
	end do
	
	call potenciali()
	call Tlakec()
	Nu=Nu+l**2*(P_real-P)*integracijski_korak/Qbar
	Eps=Eps+Nu*integracijski_korak
	l=l0*exp(Eps)
	
	scalar=exp(Nu*integracijski_korak)
	
	do i=1, N
		delci(i)%r=delci(i)%r*scalar
		delci(i)%v=delci(i)%v*scalar
	end do
	
	call potenciali()

	do i=1, N
		delci(i)%v=delci(i)%v+delci(i)%a*integracijski_korak/2.0d0
	end do
	
	call Nose_Hoover_thermostat()
	
end subroutine MD_korak

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
	call kineticna_energija()
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
		Wkin=Wkin+(delci(i)%v(1)**2+delci(i)%v(2)**2)/2
	end do
end subroutine Kineticna_energija

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
	
	P_real=vsota/(2.*l**2)+real(N)/l**2*T
end subroutine Tlakec

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

!Berendsenov barostat
subroutine Berendsen_barostat
	use parametri_simulacije
	use box
	use particles
	use meritve
	use simulacija
	implicit none
	real*8:: lambda
	real*8, parameter:: tau=0.5
	integer:: i
	
	call Tlakec()
	
	lambda=sqrt(1-integracijski_korak/tau*(P/P_real-1.))
	
	l=l*lambda
	do i=1, N
		delci(i)%r=delci(i)%r*lambda
		delci(i)%v=delci(i)%v*lambda
	end do
end subroutine Berendsen_barostat

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

!zapiše rezultate poskusa pri temperaturi T in tlaku P
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
	
	Uavg=Uavg/real((st_korakov)/korak_vzorcenja)
	U2=U2/real((st_korakov)/korak_vzorcenja)
	rhoavg=rhoavg/real((st_korakov)/korak_vzorcenja)
	Vavg=Vavg/real((st_korakov)/korak_vzorcenja)
	V2=V2/real((st_korakov)/korak_vzorcenja)
	Havg=Havg/real((st_korakov)/korak_vzorcenja)
	H2=H2/real((st_korakov)/korak_vzorcenja)
	VH=VH/real((st_korakov)/korak_vzorcenja)
	
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
	
	write(100, "(13f15.6)") T, rho_kon, s_rho_kon, U_kon, s_U_kon, H_kon, s_H_kon, cp_kon, s_cp_kon, k_kon, s_k_kon, a_kon, s_a_kon
	Uavg=0
	U2=0
	rhoavg=0
	Vavg=0
	V2=0
	Havg=0
	H2=0
	VH=0
	
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