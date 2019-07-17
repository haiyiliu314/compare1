!Module defining universal constants
!***************************************************************************
module univ_constant2
implicit none

!ee: fundamental charge 1.602176634*10^-19 C
!hbar: Plank constant 6.62607015/(2*Pi)*10^-34 J s
!alpha: fine structure constant 1/137
!rydberg: Rydberg unit of energy 13.6 eV
!a_bohr: Bohr radius 0.53 A
double precision,parameter :: Pi=3.141592653589793d0
double precision,parameter :: ee=1.602176634d0
double precision,parameter :: hbar=6.62607015d0/(2d0*Pi)
!double precision,parameter :: alpha=1d0/137.035999139d0
double precision,parameter :: rydberg=13.605693009d0
double precision,parameter :: a_bohr=0.529177210903d0
double precision,parameter :: light_c=2.99792458d0

!Effective masses in units of m_0
double precision,parameter :: mass_c=0.226733479967411d0
double precision,parameter :: mass_h=0.281802063149898d0

!reduced masses in units of m_0
double precision,parameter :: mu_ch=mass_c*mass_h/(mass_c+mass_h)

!Band edges in units of eV
double precision,parameter :: E_c=4.8894d0
double precision,parameter :: E_h=-0.1351d0

!dielectric constants
double precision,parameter :: diex_gan=8.9d0
double precision,parameter :: diez_gan=9.8d0
double precision,parameter :: die_gan=dsqrt(diex_gan*diez_gan)
double precision,parameter :: nz_aln=dsqrt(8.5d0)

!Energy units as exciton binding energy in units of eV
!length units as exciton Bohr radius in units of A
!time units as 1 ps
double precision,parameter :: E_B=rydberg*mu_ch/die_gan**2d0
double precision,parameter :: a_B=die_gan*a_bohr/mu_ch
double precision,parameter :: tau=1d0

!prefactor of SBE: -iE_{tau}=-i*tau*E_B/hbar
double precision,parameter :: E_tau_r=1000d0*ee*tau*E_B/hbar
double complex,parameter :: E_tau=dcmplx(0d0,-E_tau_r)
double complex,parameter :: E_tau2=E_tau*2d0

!Band gaps in units of E_B
double precision,parameter :: Eg_bulk=3.507d0/E_B
double precision,parameter :: E_ch=(E_c-E_h)/E_B

!Dipole constants in units of e A
!overlap of envelope functions
double precision,parameter :: fch_overlap=0.891750550800708d0 
double precision,parameter :: d_bulk=2.09d0
double precision,parameter :: d_ch=fch_overlap*d_bulk*Eg_bulk/E_ch

!External field parameters
!field strength in units of kV/cm
double precision,parameter :: Ele_source=2000d0
double precision,parameter :: w_factor=E_tau_r/E_B
double precision,parameter :: A0_factor=a_B*Ele_source*0.00001d0

!Excitation induced dephasing
!double precision,parameter :: e_correct=(4d0*rydberg*a_bohr)/(nz_aln*E_B*light_c*tau)*(d_ch/a_B)**2d0

!Ratio of E_r/E_i
double precision,parameter :: ratio_et=(0.4d0*d_ch*rydberg*a_bohr)/(nz_aln*Ele_source*light_c*tau*(a_B**2d0))

!Detuning and dephasing
double precision,parameter :: gamma_1=0.000d0/E_B
double precision,parameter :: gamma_2=0.005d0/E_B
double precision,parameter :: gamma_qm=5d0*gamma_2
double precision,parameter :: E_shift=0d0/E_B
double complex,parameter :: f_detune=dcmplx(0d0,-gamma_1)
double complex,parameter :: c_detune=dcmplx(0d0,-gamma_qm)
double complex,parameter :: g_detune=dcmplx(0d0,gamma_2)

!units of dipole coupling
double precision,parameter :: E_dot_dch=d_ch*Ele_source*0.00001d0/E_B

end module univ_constant2

!********************************************************************

!Main program
!********************************************************************
program gan_sbe_3

use univ_constant2
implicit none


!Field correction
double precision :: field_correct,ratio_correct
!********************************************************************

!********************************************************************
!Grid y_j(:)
!Range of |k|: 0-y_max in units of 1/a_B
!Interval: delta_y
!Number of points: n_y
!Finer grid for diagonal elements of Coulomb potential: m_y
double precision,pointer :: y_j(:)
double precision :: y_max,delta_y
integer :: n_y,m_y

!Intermediate vectors from y-grid
!yj_sq=y^2
!A0_y=y*A0
double precision,pointer :: yj_sq(:),A0_y(:)
!Ech_y2=E_ch+y^2
double complex,pointer :: Ech_y2(:),Eqm_y2(:)
!********************************************************************

!********************************************************************
!Coulomb matrices
!Factor for Coulomb interaction
!g_kpara as a function of the sample points y_para
!number of points: n_sample
double precision,pointer :: y_para(:)
double precision,pointer :: gch_para(:)
double precision,pointer :: gcc_para(:),ghh_para(:)
integer :: n_sample

!order_polar: cut-off order of polar expansion
!num_ml: the number of monolayer
!str_ml: num_ml to string
!char_head,char_tail: Folder name for g_k factor
integer :: order_polar
integer :: num_ml
character(1) ::str_ml
character(30),parameter :: char_head='/qstlhome/liuhai/run/Add_confinement/code/'
character(13),parameter :: char_taily='ml/y_para.txt'
character(4),parameter :: char_midg='ml/g'
character(9),parameter :: char_tailg='_para.txt'
character(2),parameter :: char_ch='ch'
character(2),parameter :: char_cc='cc'
character(2),parameter :: char_hh='hh'
character(44) :: yfile_name
character(46) :: gfile_name
character(len=100)  ::format_V,format_V1
character(80)  ::list_file
double precision,pointer :: Vcoul_ch(:,:,:)
double precision,pointer :: Vcoul_cc(:,:,:),Vcoul_hh(:,:,:)
!********************************************************************


!********************************************************************
!CPU clock timing
double precision :: tic_cpu,toc_cpu
!********************************************************************


!********************************************************************
!Runge Kutta
double complex,pointer :: p_ch1(:,:),p_ch0(:,:),p_chm(:,:)
double complex,pointer :: g_qm1(:,:),g_qm0(:,:),g_qmm(:,:)
double complex,pointer :: f_e1(:,:),f_e0(:,:),f_em(:,:)

double complex,pointer :: kp_ch(:,:),kg_qm(:,:)
double complex,pointer :: kf_e(:,:)

double complex,pointer :: omega_ch(:,:),omega_qm(:,:)
double complex,pointer :: omega_cnh(:,:)

double precision :: E_t,At_sq,a_shape,a_int
double complex :: At_plus,At_minus,Atp_sq,Atm_sq

double precision :: delta_t,hdelta_t,thirdelta_t,sxdelta_t
double precision :: f_max,w_broad
double precision :: t0,thf,t1
integer :: m_max,p_order,order_p2
integer :: time_step
integer :: j_rk4, i1
!********************************************************************

!********************************************************************
!Sweeping the frequency
!field frequency as photon energy in units of eV
double precision :: Efr_source0,Efr_source,w_source
double precision :: sweep_w,sweep_step,sweep_end
!units of vector potential
double precision :: A_0,A0_sq
!********************************************************************
integer  ::readin

!********************************************************************
!Get the y-grid
y_max=30d0
n_y=200
m_y=50
allocate(y_j(n_y))
call half_inte_gri(y_j,delta_y,y_max,n_y)

!calculate the kinetic energy and light-matter coupling
allocate(yj_sq(n_y),A0_y(n_y))
yj_sq(:)=y_j(:)**2d0

allocate(Ech_y2(n_y),Eqm_y2(n_y))
Ech_y2(:)=(E_ch+E_shift)+yj_sq(:)
Eqm_y2(:)=Ech_y2(:)+c_detune

!********************************************************************

!Calculate the Coulomb matrices
!Load the sample values for the factor g(q)
!**********************************************************************
n_sample=101
allocate(gch_para(n_sample))
allocate(gcc_para(n_sample),ghh_para(n_sample))
allocate(y_para(n_sample))
num_ml=2
write(str_ml,'(i1)'),num_ml

yfile_name='y_para.txt'
open(unit=11,file=yfile_name)
read(11,*) y_para
close(11)

gfile_name='gch_para.txt'
open(unit=11,file=gfile_name)
read(11,*) gch_para
close(11)

gfile_name='gcc_para.txt'
open(unit=11,file=gfile_name)
read(11,*) gcc_para
close(11)

gfile_name='ghh_para.txt'
open(unit=11,file=gfile_name)
read(11,*) ghh_para
close(11)

!*********************************************************************

order_polar=2 !cutoff at order_polar-1
allocate(Vcoul_ch(n_y,n_y,order_polar))
allocate(Vcoul_cc(n_y,n_y,order_polar),Vcoul_hh(n_y,n_y,order_polar))

call coulomb_mat(Vcoul_ch,Vcoul_cc,Vcoul_hh,order_polar, &
	& delta_y,y_j,n_y,m_y, &
	& y_para,gch_para,gcc_para,ghh_para,n_sample)
  write(list_file, '(A)') 'Vcoul_ch.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', 1, '(", ",SE24.16e3))'
  write(format_V1, '(A12, I4, A18)')   '(SE24.16e3, ', n_y, '(", ",SE24.16e3))'
  do i1 = 1, n_y
      write(700, format_V) real(Vcoul_ch(:,i1,1))
  end do
Vcoul_cc(:,:,:)=(Vcoul_cc(:,:,:)+Vcoul_hh(:,:,:))/2d0
!===================================================================================================================================================================================!                 
! Vcoul_cc(:,:,:) = 0.0d0
! Vcoul_ch(:,:,:) = 0.0d0
!===================================================================================================================================================================================
open(unit=20,file='univ_constants.txt')
	write(20,*),'Effective masses:'
	write(20,*),'m_c=',mass_c,' m_0'
	write(20,*),'m_h=',mass_h,' m_0'

	write(20,*),'Reduced masses:'
	write(20,*),'mu_ch=',mu_ch,' m_0'

	write(20,*),'Band edges:'
	write(20,*),'E_c=',E_c,' eV'
	write(20,*),'E_h=',E_h,' eV'
	
	write(20,*),'Dielectric constants:'
	write(20,*),'epsilon_x=',diex_gan
	write(20,*),'epsilon_z=',diez_gan
	write(20,*),'epsilon_ave=',die_gan

	write(20,*),'Units:'
	write(20,*),'E_B = ',E_B,' eV'	
	write(20,*),'a_B = ',a_B,' A'	
	write(20,*),'tau = ',tau,' ps'
	write(20,*),'E_tau = tau*E_B/hbar=',E_tau_r
	
	write(20,*),'Eg_bulk = ',Eg_bulk,'E_B'
	write(20,*),'E_ch = ',E_ch,'E_B'

	write(20,*),'d_bulk = ',d_bulk,'e A'
	write(20,*),'d_ch = ',d_ch,'e A'

	write(20,*),'Ele_source = ',Ele_source,'kV/cm'

	write(20,*),'gamma_2 = ',gamma_2,'E_B'
	write(20,*),'gamma_qm = ',gamma_qm,'E_B'

	write(20,*),'m_max=',order_polar-1
	write(20,*),'y_max=',y_max
	write(20,*),'y-grid: ',n_y
	write(20,*),'fine grid: ',m_y

!*********************************************************************
!Determine the timestep
!f_max=E_tau_r*(E_cl+E_shift+r_hl*y_max**2d0)    !Edit_label
!delta_t=1d0/f_max

w_broad=gamma_2/(2d0*Pi) !Broadening of the driving field in frequency domain in units of E_B
!w_broad=10d0 
f_max=E_tau_r*w_broad
write(*,*) f_max
!delta_t=1/(E_tau_r*y_max**2d0)
!delta_t=1/(800d0*f_max) 
delta_t=Pi/(4000d0*E_tau_r) 

hdelta_t=delta_t/2d0
thirdelta_t=delta_t/3d0
sxdelta_t=delta_t/6d0

time_step=int(6d0/(gamma_2*E_tau_r*delta_t))           !Errors depend on the time_step

if(gamma_2>w_broad) then
	time_step=int(6d0/(w_broad*E_tau_r*delta_t))           !Errors depend on the time_step
end if
!time_step = 10
!*********************************************************************



	order_p2=2*order_polar
	p_order=order_p2-1
	m_max=order_polar-1

	allocate(p_ch1(n_y,p_order),p_ch0(n_y,p_order),p_chm(n_y,p_order))
	allocate(g_qm1(n_y,p_order),g_qm0(n_y,p_order),g_qmm(n_y,p_order))
	
	allocate(f_e1(n_y,p_order),f_e0(n_y,p_order),f_em(n_y,p_order))

	allocate(kp_ch(n_y,p_order),kg_qm(n_y,p_order))
	allocate(kf_e(n_y,p_order))

	allocate(omega_ch(n_y,p_order),omega_qm(n_y,p_order))
	allocate(omega_cnh(n_y,p_order))
  write(list_file, '(A)') 'omega_sweep.dat'
  open(unit=11,file=list_file)
  write(list_file, '(A)') 'density.dat'
  open(unit=14,file=list_file)
  write(list_file, '(A)') 'test.dat'
  open(unit=15,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', 1, '(", ",SE24.16e3))'
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', n_y, '(", ",SE24.16e3))'
!=============================================================================================================================================================

open(unit = 100, file = 'fort.42', status = 'old', action = 'read')
read(100,*) readin

sweep_end=25d-3
sweep_step=0.6d-3
Efr_source0=5.018528494722949d0/2d0
sweep_w=-30d-3+dble(readin)*sweep_step
ratio_correct=ratio_et*delta_y

call cpu_time(tic_cpu)



    Efr_source=Efr_source0+sweep_w

    w_source=Efr_source*w_factor
    a_int=w_source*hdelta_t
    A_0=A0_factor/Efr_source
    A0_y(:)=y_j(:)*A_0

    A0_sq=A_0**2d0

	p_ch0(:,:)=0d0
	g_qm0(:,:)=0d0
	f_e0(:,:)=0d0

	omega_qm(:,:)=0d0
	omega_ch(:,:)=0d0
	omega_cnh(:,:)=0d0

	t0=0d0
	a_shape=0d0
    	field_correct=0d0
	call field_t(E_t,At_sq,At_plus,Atp_sq,f_max,t0,a_shape,a_int,field_correct,w_source)

	do j_rk4=1,time_step

		!k1=F(z(t),t)
		call rhs_sbe(kp_ch,kg_qm,kf_e,p_ch0,g_qm0,f_e0,order_polar,m_max,p_order,order_p2, &
		     & omega_ch,omega_qm,omega_cnh,n_y, &
                     & Ech_y2,Eqm_y2,A0_y,E_t,A0_sq,At_sq,At_plus,Atp_sq)

		p_ch1(:,:)=p_ch0(:,:)+sxdelta_t*kp_ch(:,:)
		g_qm1(:,:)=g_qm0(:,:)+sxdelta_t*kg_qm(:,:)
		f_e1(:,:)=f_e0(:,:)+sxdelta_t*kf_e(:,:)

		thf=t0+hdelta_t

        	!field_correct=ratio_correct*real(sum(kp_ch(:,order_polar)*y_j(:)))

		call field_t(E_t,At_sq,At_plus,Atp_sq,f_max,thf,a_shape,a_int,field_correct,w_source)

		p_chm(:,:)=p_ch0(:,:)+hdelta_t*kp_ch(:,:)
		g_qmm(:,:)=g_qm0(:,:)+hdelta_t*kg_qm(:,:)
		f_em(:,:)=f_e0(:,:)+hdelta_t*kf_e(:,:)

		!k2=F(z(t)+hdelta_t*k1,t+hdelta_t)

		call coulomb_sum(p_chm,g_qmm,f_em,order_polar,p_order,n_y, &   
			& Vcoul_ch,Vcoul_cc, &
			& omega_ch,omega_qm,omega_cnh)

		call rhs_sbe(kp_ch,kg_qm,kf_e,p_chm,g_qmm,f_em,order_polar,m_max,p_order,order_p2, &
		     & omega_ch,omega_qm,omega_cnh,n_y, &
                     & Ech_y2,Eqm_y2,A0_y,E_t,A0_sq,At_sq,At_plus,Atp_sq)

		p_ch1(:,:)=p_ch1(:,:)+thirdelta_t*kp_ch(:,:)
		g_qm1(:,:)=g_qm1(:,:)+thirdelta_t*kg_qm(:,:)
		f_e1(:,:)=f_e1(:,:)+thirdelta_t*kf_e(:,:)

		p_chm(:,:)=p_ch0(:,:)+hdelta_t*kp_ch(:,:)
		g_qmm(:,:)=g_qm0(:,:)+hdelta_t*kg_qm(:,:)
		f_em(:,:)=f_e0(:,:)+hdelta_t*kf_e(:,:)

		!k3=F(z(t)+hdelta_t*k2,t+hdelta_t)

		call coulomb_sum(p_chm,g_qmm,f_em,order_polar,p_order,n_y, &   
			& Vcoul_ch,Vcoul_cc, &
			& omega_ch,omega_qm,omega_cnh)

		call rhs_sbe(kp_ch,kg_qm,kf_e,p_chm,g_qmm,f_em,order_polar,m_max,p_order,order_p2, &
		     & omega_ch,omega_qm,omega_cnh,n_y, &
                     & Ech_y2,Eqm_y2,A0_y,E_t,A0_sq,At_sq,At_plus,Atp_sq)
			
		p_ch1(:,:)=p_ch1(:,:)+thirdelta_t*kp_ch(:,:)
		g_qm1(:,:)=g_qm1(:,:)+thirdelta_t*kg_qm(:,:)
		f_e1(:,:)=f_e1(:,:)+thirdelta_t*kf_e(:,:)

		t1=t0+delta_t

		call field_t(E_t,At_sq,At_plus,Atp_sq,f_max,t1,a_shape,a_int,field_correct,w_source)

		p_chm(:,:)=p_ch0(:,:)+delta_t*kp_ch(:,:)
		g_qmm(:,:)=g_qm0(:,:)+delta_t*kg_qm(:,:)
		f_em(:,:)=f_e0(:,:)+delta_t*kf_e(:,:)

		!k4=F(z(t)+delta_t*k3,t+delta_t)

		call coulomb_sum(p_chm,g_qmm,f_em,order_polar,p_order,n_y, &   
			& Vcoul_ch,Vcoul_cc, &
			& omega_ch,omega_qm,omega_cnh)

		call rhs_sbe(kp_ch,kg_qm,kf_e,p_chm,g_qmm,f_em,order_polar,m_max,p_order,order_p2, &
		     & omega_ch,omega_qm,omega_cnh,n_y, &
                     & Ech_y2,Eqm_y2,A0_y,E_t,A0_sq,At_sq,At_plus,Atp_sq)
		
		p_ch1(:,:)=p_ch1(:,:)+sxdelta_t*kp_ch(:,:)
		g_qm1(:,:)=g_qm1(:,:)+sxdelta_t*kg_qm(:,:)
		f_e1(:,:)=f_e1(:,:)+sxdelta_t*kf_e(:,:)

		t0=t1
		p_ch0(:,:)=p_ch1(:,:)
		g_qm0(:,:)=g_qm1(:,:)
		f_e0(:,:)=f_e1(:,:)

		call coulomb_sum(p_ch0,g_qm0,f_e0,order_polar,p_order,n_y, &
			& Vcoul_ch,Vcoul_cc, &
			& omega_ch,omega_qm,omega_cnh)


    write(15,format_V)   E_dot_dch*E_t
	end do
    write(14,format_V) dble(sum(f_e0(:,order_polar)*y_j(:)))
!    write(*,*) sweep_w
    write(11,format_V) sweep_w



 close(11)
 close(14)
 close(15)
call cpu_time(toc_cpu)

write(20,*),'Total time:',toc_cpu-tic_cpu,'s'
write(20,*),'Number of Steps: ',time_step
write(20,*),'Time for each step:',(toc_cpu-tic_cpu)/time_step,'s'
write(20,*),'dt in SBE:',delta_t, 'ps'
write(20,*),'Total evolution time:',delta_t*time_step,'ps'
write(20,*),'Pulse half width',0.8326d0/f_max,'ps'

close(20)

contains

subroutine coulomb_sum(p_ch,g_qm,f_e,order_polar,p_order,n_y, &
			& Vcoul_ch,Vcoul_cc, &
			& omega_ch,omega_qm,omega_cnh)

	implicit none
	!INPUT
	integer :: order_polar,p_order,n_y
	double complex,pointer :: p_ch(:,:),g_qm(:,:)
	double complex,pointer :: f_e(:,:)

	double precision,pointer :: Vcoul_ch(:,:,:)
	double precision,pointer :: Vcoul_cc(:,:,:)

	!OUTPUT
	double complex,pointer :: omega_ch(:,:),omega_qm(:,:)
	double complex,pointer :: omega_cnh(:,:)

	integer :: j,j_coul

	do j=1,p_order

		j_coul=abs(j-order_polar)+1

		omega_ch(:,j)=matmul(Vcoul_ch(:,:,j_coul),p_ch(:,j))

		omega_qm(:,j)=matmul(Vcoul_ch(:,:,j_coul),g_qm(:,j))

		omega_cnh(:,j)=matmul(Vcoul_cc(:,:,j_coul),f_e(:,j))

	end do

end subroutine coulomb_sum

subroutine field_t(E_t,At_sq,At_plus,Atp_sq,f_max,t,a_shape,a_int,field_correct,w_source)

	implicit none

	!Calculate the dipole-coupling term & terms associated with the vector potential
	
	!OUTPUT
	double precision :: E_t,At_sq
	double complex :: At_plus,Atp_sq

	!INPUT
	double precision :: t
    	double precision :: field_correct,w_source

	!prefactors
	double precision :: f_max,a_int
	
	!intermediate variables
	double precision :: a_shape

    	E_t=dexp(-(f_max*t-3d0)**2d0)*dcos(w_source*t) !-field_correct	!Test_label

	a_shape=a_shape-E_t*a_int
	!a_shape=-dsin(w_drive*t)	!Test_label

	At_sq=a_shape**2d0
	
	At_plus=dcmplx(1d0,0d0)*a_shape
	Atp_sq=At_plus**2d0

end subroutine field_t

subroutine rhs_sbe(kp_ch,kg_qm,kf_e,p_ch,g_qm,f_e,order_polar,m_max,p_order,order_p2, &
		     & omega_ch,omega_qm,omega_cnh,n_y, &
                     & Ech_y2,Eqm_y2,A0_y,E_t,A0_sq,At_sq,At_plus,Atp_sq)
		
	implicit none

	!Calulate the right-hand side of SBE

	!OUTPUT
	double complex,pointer :: kp_ch(:,:),kg_qm(:,:)
	double complex,pointer :: kf_e(:,:)

	!INPUT
	double complex,pointer :: p_ch(:,:),g_qm(:,:)
	double complex,pointer :: f_e(:,:)

	double complex,pointer :: omega_ch(:,:),omega_qm(:,:)
	double complex,pointer :: omega_cnh(:,:)

	double complex,pointer :: Ech_y2(:),Eqm_y2(:)
	double precision,pointer :: A0_y(:)	

	double precision :: E_t,A0_sq,At_sq
	double complex :: At_plus,Atp_sq

	integer :: order_polar,m_max
	integer :: p_order,order_p2
	integer :: n_y

	!INTERMEDIATE VARIABLES
	double complex,dimension(n_y,p_order) :: comega_ch

	double complex,dimension(n_y,p_order) :: domega_ch
	double complex,dimension(n_y,p_order) :: media_qm1,media_qm2
	double complex,dimension(n_y,p_order) :: media_ch1,media_ch2
	double complex,dimension(n_y,p_order) :: media_e


    	double complex,dimension(n_y) :: mvec_ch
   	double complex,dimension(n_y) :: mvec_chkap,mvec_chkam

	!Variable for do loop
	integer :: j,m_pr,mm_pr
	integer :: j_med1

	domega_ch(:,:)=omega_ch(:,:)

	domega_ch(:,order_polar)=domega_ch(:,order_polar)+E_dot_dch*E_t   

	comega_ch(:,:)=conjg(domega_ch(:,p_order:1:-1))

	media_qm1(:,:)=0d0
	media_qm2(:,:)=0d0

	media_ch1(:,:)=0d0
	media_ch2(:,:)=0d0

	media_e(:,:)=0d0

	do j=1,p_order
		
		do m_pr=max(1,j-m_max),min(p_order,j+m_max)

			mm_pr=j-m_pr+order_polar

			media_qm1(:,j)=media_qm1(:,j)+omega_cnh(:,m_pr)*g_qm(:,mm_pr)
			media_qm2(:,j)=media_qm2(:,j)+omega_qm(:,m_pr)*f_e(:,mm_pr)

			media_ch1(:,j)=media_ch1(:,j)+omega_cnh(:,m_pr)*p_ch(:,mm_pr)
			media_ch2(:,j)=media_ch2(:,j)+domega_ch(:,m_pr)*f_e(:,mm_pr)

			media_e(:,j)=media_e(:,j)+p_ch(:,m_pr)*comega_ch(:,mm_pr)

		end do

	end do


	media_e(:,:)=media_e(:,:)-conjg(media_e(:,p_order:1:-1))

    	mvec_ch(:)=Ech_y2(:)+A0_sq*At_sq

    	mvec_chkap(:)=A0_y(:)*At_plus
    	mvec_chkam(:)=conjg(mvec_chkap(:))
 !       mvec_chkam(:) = 0.0d0

	do j=1,p_order

		kp_ch(:,j)=mvec_ch(:)*p_ch(:,j)
		kg_qm(:,j)=Eqm_y2(:)*g_qm(:,j)		
                         
		j_med1=j-1
		if(j_med1>=1 .AND. j_med1<=p_order) then
			kp_ch(:,j)=kp_ch(:,j)+mvec_chkam(:)*p_ch(:,j_med1)
			!kg_qm(:,j)=kg_qm(:,j)+mvec_chkam(:)*g_qm(:,j_med1)
		end if
		
		j_med1=j+1		
		if(j_med1>=1 .AND. j_med1<=p_order) then
			kp_ch(:,j)=kp_ch(:,j)+mvec_chkap(:)*p_ch(:,j_med1)
			!kg_qm(:,j)=kg_qm(:,j)+mvec_chkap(:)*g_qm(:,j_med1)
		end if

	end do


    	kp_ch(:,:)=E_tau*(kp_ch(:,:)+c_detune*g_qm(:,:)-domega_ch(:,:)-media_ch1(:,:)+media_ch2(:,:))

    	kg_qm(:,:)=E_tau*(kg_qm(:,:)+g_detune*p_ch(:,:)-omega_qm(:,:)-media_qm1(:,:)+media_qm2(:,:))
  
	kf_e(:,:)=E_tau2*media_e(:,:)  !Delete dephasing


end subroutine rhs_sbe

subroutine half_inte_gri(y,delta_y,y_max,n)

	implicit none

	!Grid starting at y=0

	!OUTPUT
	double precision,pointer :: y(:)
	double precision :: delta_y

	!INPUT
	double precision :: y_max
	integer :: n

	!Variable for do loop
	integer :: j

	delta_y=y_max/n

	do j=1,n
		y(j)=delta_y*(j-0.5d0)
	end do

end subroutine half_inte_gri

subroutine coulomb_mat(Vcoul_ch,Vcoul_cc,Vcoul_hh,order_polar, &
	& delta_y,y_j,n_y,m_y, &
	& y_para,gch_para,gcc_para,ghh_para,n_sample)

	implicit none

	!Calculate the Coulomb matrices

	!OUTPUT
	double precision,pointer :: Vcoul_ch(:,:,:)
	double precision,pointer :: Vcoul_cc(:,:,:),Vcoul_hh(:,:,:)

	!INPUT
	!y-grid
	double precision,pointer :: y_j(:)
	double precision :: delta_y
	integer :: n_y
	!cut-off order of polar expansion
	integer :: order_polar

	!Grid of modified Gauss-Chebyshev quadratures
	double precision,pointer :: theta(:),w_cheby(:)
	double precision,pointer :: fch_cheby(:,:)
	double precision,pointer :: fcc_cheby(:,:),fhh_cheby(:,:)
	integer :: n_cheby

	!Finer grid for diagonal element of Coulomb potential
	double precision,pointer :: v_l(:)
	double precision,pointer :: vch_small(:,:)
	double precision,pointer :: vcc_small(:,:),vhh_small(:,:)
	integer :: m_y

	!Intermediate variables
	double precision :: y_shif,m_pi
	integer :: circ_m,j1,j2		
	
	!Factor for Coulomb interaction
	!g_kpara as a function of the sample points y_para
	!number of points: n_sample
	double precision,pointer :: y_para(:)
	double precision,pointer :: gch_para(:)
	double precision,pointer :: gcc_para(:),ghh_para(:)
	integer :: n_sample

	!Get the grid of modified Gauss-Chebyshev quadratures
	!f_cheby to be the function values from the Coulomb potential
	n_cheby=100		
	allocate(theta(n_cheby),w_cheby(n_cheby))
	call modified_chebyshev_grid(theta,w_cheby,n_cheby)	
	

	!Initialization of the matrices
	Vcoul_ch(:,:,:)=0d0
	Vcoul_cc(:,:,:)=0d0
	Vcoul_hh(:,:,:)=0d0
	
	allocate(fch_cheby(n_cheby,order_polar))
        allocate(fcc_cheby(n_cheby,order_polar),fhh_cheby(n_cheby,order_polar))


	m_pi=delta_y/Pi    !Edit_label

	!nondiagonal matrix elements

	do j1=1,n_y
		do j2=j1+1,n_y
			call coulomb_vyy(fch_cheby,fcc_cheby,fhh_cheby,order_polar, &
 & y_j(j2),y_j(j1),theta,n_cheby,y_para,gch_para,gcc_para,ghh_para,n_sample)

			do circ_m=1,order_polar
				Vcoul_ch(j2,j1,circ_m)=m_pi*sum(w_cheby(:)*fch_cheby(:,circ_m))  !Coulomb matrix is symmetric
				Vcoul_ch(j1,j2,circ_m)=y_j(j2)*Vcoul_ch(j2,j1,circ_m)
				Vcoul_ch(j2,j1,circ_m)=y_j(j1)*Vcoul_ch(j2,j1,circ_m)

				Vcoul_cc(j2,j1,circ_m)=m_pi*sum(w_cheby(:)*fcc_cheby(:,circ_m))  !Coulomb matrix is symmetric
				Vcoul_cc(j1,j2,circ_m)=y_j(j2)*Vcoul_cc(j2,j1,circ_m)
				Vcoul_cc(j2,j1,circ_m)=y_j(j1)*Vcoul_cc(j2,j1,circ_m)

				Vcoul_hh(j2,j1,circ_m)=m_pi*sum(w_cheby(:)*fhh_cheby(:,circ_m))  !Coulomb matrix is symmetric
				Vcoul_hh(j1,j2,circ_m)=y_j(j2)*Vcoul_hh(j2,j1,circ_m)
				Vcoul_hh(j2,j1,circ_m)=y_j(j1)*Vcoul_hh(j2,j1,circ_m)

			end do
		end do
	end do
	

	!digonal matrix elements
	allocate(v_l(m_y))
	allocate(vch_small(m_y,order_polar))
	allocate(vcc_small(m_y,order_polar),vhh_small(m_y,order_polar))

	y_shif=(0.5d0+0.5d0/m_y)*delta_y

	m_pi=m_pi/m_y

	do j1=1,n_y
		!Get the grid for interval labeled by j1
		call half_integer_grid(v_l,y_j(j1)-y_shif,delta_y,m_y)

		do j2=1,m_y
			call coulomb_vyy(fch_cheby,fcc_cheby,fhh_cheby,order_polar, &
 & v_l(j2),y_j(j1),theta,n_cheby,y_para,gch_para,gcc_para,ghh_para,n_sample)

			do circ_m=1,order_polar
				vch_small(j2,circ_m)=v_l(j2)*sum(w_cheby(:)*fch_cheby(:,circ_m))
				vcc_small(j2,circ_m)=v_l(j2)*sum(w_cheby(:)*fcc_cheby(:,circ_m))
				vhh_small(j2,circ_m)=v_l(j2)*sum(w_cheby(:)*fhh_cheby(:,circ_m))
			end do
		end do

!write(*,*) w_cheby
		do circ_m=1,order_polar
			Vcoul_ch(j1,j1,circ_m)=sum(vch_small(:,circ_m))*m_pi
			Vcoul_cc(j1,j1,circ_m)=sum(vcc_small(:,circ_m))*m_pi
			Vcoul_hh(j1,j1,circ_m)=sum(vhh_small(:,circ_m))*m_pi
		end do
	end do	

end subroutine coulomb_mat


subroutine modified_chebyshev_grid(theta,w,n)

	implicit none

	!Grid of modified Gauss-Chebyshev quadratures

	!OUTPUT
	double precision,pointer :: theta(:),w(:)

	!INPUT	
	integer :: n

	!Intermediate variables
	double precision :: delta_x,dx_j
	integer :: j

	delta_x=2d0*Pi/(n+1d0)

	do j=1,n
		dx_j=delta_x*j
		theta(j)=(dx_j-dsin(dx_j))/2d0
		w(j)=delta_x*(1d0-dcos(dx_j))
	end do

end subroutine modified_chebyshev_grid

subroutine coulomb_vyy(fch_cheby,fcc_cheby,fhh_cheby,order_polar, &
 & y1,y2,theta,n,y_para,gch_para,gcc_para,ghh_para,n_sample)

	implicit none

	!Calculate the Coulom potential

	!OUTPUT
	double precision,pointer :: fch_cheby(:,:)
	double precision,pointer :: fcc_cheby(:,:),fhh_cheby(:,:)

	!INPUT
	double precision,pointer :: theta(:)
	double precision :: y1,y2
	integer :: n,order_polar,circ_m

	!Factor for Coulomb interaction
	!g_kpara as a function of the sample points y_para
	!number of points: n_sample
	double precision,pointer :: y_para(:)
	double precision,pointer :: gch_para(:)
	double precision,pointer :: gcc_para(:),ghh_para(:)
	integer :: n_sample

	!Intermediate variables
	double precision :: y_sq,y12,y
	double precision :: gch_y,gcc_y,ghh_y

	!Variable for do loop
	integer :: j

	y_sq=y1*y1+y2*y2
	y12=2d0*y1*y2

	do j=1,n
		y=dsqrt(y_sq-y12*dcos(theta(j)))
		call interp_gy(y_para,gch_para,gcc_para,ghh_para,n_sample,y,gch_y,gcc_y,ghh_y)

		do circ_m=0,order_polar-1
			fch_cheby(j,circ_m+1)=gch_y*dcos(circ_m*theta(j))/y
			fcc_cheby(j,circ_m+1)=gcc_y*dcos(circ_m*theta(j))/y
			fhh_cheby(j,circ_m+1)=ghh_y*dcos(circ_m*theta(j))/y
		end do
	end do
!write(*,*) y
end subroutine coulomb_vyy

subroutine interp_gy(y_para,gch_para,gcc_para,ghh_para,n_sample,y,gch_y,gcc_y,ghh_y)

	implicit none

	!Get the Factor g(y) for Coulomb interaction by interpolation

	!OUTPUT
	double precision :: gch_y,gcc_y,ghh_y

	!INPUT
	double precision :: y
	double precision,pointer :: y_para(:)
	double precision,pointer :: gch_para(:)
	double precision,pointer :: gcc_para(:),ghh_para(:)
	integer :: n_sample

	
	!Variable indices for interpolation
	integer :: ind_r,ind_l
	double precision :: y_ratio

	ind_r=1	
	do while(y_para(ind_r)<y)
		ind_r=ind_r+1
	end do
	ind_l=ind_r-1

	y_ratio=(y-y_para(ind_l))/(y_para(ind_r)-y_para(ind_l))

	gch_y=gch_para(ind_l)+y_ratio*(gch_para(ind_r)-gch_para(ind_l))
	gcc_y=gcc_para(ind_l)+y_ratio*(gcc_para(ind_r)-gcc_para(ind_l))
	ghh_y=ghh_para(ind_l)+y_ratio*(ghh_para(ind_r)-ghh_para(ind_l))

end subroutine interp_gy

subroutine half_integer_grid(v,v0,delta_y,m)

	implicit none

	!Grid starting at v0

	!OUTPUT
	double precision,pointer :: v(:)
	double precision :: delta_v

	!INPUT
	double precision :: delta_y
	double precision :: v0
	integer :: m

	!Variable for do loop
	integer :: j

	delta_v=delta_y/m

	do j=1,m
		v(j)=v0+delta_v*j		
	end do

end subroutine half_integer_grid

end program gan_sbe_3
!********************************************************************
