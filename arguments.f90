  module constants  
  !07/03/2017 Creation
  !Initialization and calculation of all the constants and global variables 
  use N_of_grids
  implicit none
          
  !t_end: whole time
!  integer                                      ::m(2*Nm_o+1) = 0
!  !m: order of circular harmonics

  !ymax: maximum value for y   dy: length of one step for y grid

!                                                 sigmat_A = 0.2d0*scale1, tstart_A = -0.5d0*scale1  
  !sigmat: sigma for gaussian pulse(E(t) = exp(-t^2/sigmat^2))(ps)   tstart: starting point for
  !time (ps) 
  double precision, parameter                  ::pi = 4.0d0*datan(1.0d0)        
  !physical parameters
  double precision, parameter                  ::hbar = 0.6582119514d0, e = 1.60217662d-19, f_s =  0.0072973525693d0, eps_r = 9.33916484488844d0, c_vac = 299792458.0d-12
  !physical parameters  unit of hbar: meV * ps, f_s: fine-structure constant, eps_r: relative permittivity, c_vac: speed of light
  double precision, parameter                  ::E_B = 19.599414370216266d0, Eg = 4.8894d3+0.1351d3, a_B = eps_r*0.52917721067d0/0.125643061344944d0!a_B = hbar*c_vac*f_s/(2.0d0*eps_r*E_B)
  !E_B: binding energy(meV)   gamma: dephasing factor (meV) Eg: ground state energy
  !decay_m: dephasing caused by quantum memory unit: meV
  double precision, parameter                  ::omega_1s = (0d0)/hbar   
                                                 !A_freq_para = (Eg-4d0*E_B+2d0*E_B+12d0*E_B*(20d0/50d0 - 0.5d0))/2d0/hbar
                                                 
!                                                 A_freq_para = omega_1s
!                                                dipole: unit: m*e(elementary charge)
         
  !length of one time step

!-------------------LAPACK-----------------------
  integer, parameter                           ::LDA = Ny, LDVL = Ny, LDVR = Ny
  integer, parameter                           ::LWMAX = 10000
  integer                                      ::INFO, LWORK, LWORK1
  double precision                             ::RWORK(2*Ny)
  complex*16                                   ::A( LDA, Ny ), VL( LDVL, Ny ), VR( LDVR, Ny ),&
                                                 W( Ny ), WORK( LWMAX ), VL1( LDVL, Ny ), &
                                                 VR1( LDVR, Ny ), W1( Ny ), A1( LDA, Ny ), &
                                                 VR_temp(LDVR, Ny), VL_temp(LDVL, Ny), &
                                                 VL2( LDVL, Ny ), VR2( LDVR, Ny ), &
                                                 W2( Ny ), A2( LDA, Ny ), &
                                                 VL3( LDVL, Ny ), VR3( LDVR, Ny ), &
                                                 W3( Ny ), A3( LDA, Ny )
  character                                    ::YES = 'V',NO = 'N'
!-------------------end LAPACK-----------------------
  double precision                             ::sigmat = 0.7d0, tstart = -3.0d0,&
                                                 sigmat_A = 0.05d0
  double precision                             ::t_end = 0.0d0
  double precision                             ::ymax = 10.0d0, dy = 0d0, gamma = 60.0d0/E_B, dy1(Ny), dy_sec(N_sec)
  double precision                             ::dt = 0d0, decay_m =0.0d0/E_B
  double precision                             ::tstart_A = 0d0 
  double precision                             ::E_excit = 1.0d-3, shift = (Eg - omega_1s*hbar)/E_B,A_excit = 0.0d0
  !E_excit: excitation level(unit: binding energy)  shift: shift caused by rotation frame 
  !(unit: binding energy)
  !A_excit: unit: V*ps/m
  double precision                             ::coul_mat(Nm_o+1, Ny, Ny) = 0.0d0, &
                                                 coul_mat_cc(Nm_o+1, Ny, Ny) = 0.0d0, &
                                                 coul_mat_ch(Nm_o+1, Ny, Ny) = 0.0d0, &
                                                 coul_mat_hh(Nm_o+1, Ny, Ny) = 0.0d0, &
                                                 coul_mat_f(Nm_o+1, Ny, Ny) = 0.0d0, &
                                                 Etemp(Ny), Etemp1(Ny), Emax, dipole(Ny) = 1.30086607416965d0, y_sec(N_sec) !(A)
  !coul_mat: Coulomb matrix(non-symmetric)(unit: binding energy) in 2D Et: electrical field 
  !for excitation
  double precision                             ::y(Ny)=0.0d0, delta_1s = 0.0d0
  !y: grid for y
  double precision                             ::y_fine(N_fine) = 0.0d0, dy_fine, dy_fine_sec(N_sec), A_freq_para = 0.0d0
  !y_fine: finer grid for removal of singularity   dy_fine: length of one step of finer grid for
  !removal of ringularity  f:density
  integer                                      ::Ny1_help(N_sec+1) = 0
  complex*16                                   ::p(2*Nm_o+1, Ny) = 0.0d0, ii = (0.0d0, 1.0d0), &
                                                 f(2*Nm_o+1, Ny) = 0.0d0, coup(Ny), &
                                                 decay(2*Nm_o+1, Ny) = 0.0d0, &
                                                 p_proj(Ny), p_proj1(Ny), &
                                                 p_proj2(Ny), p_proj11(Ny), &
                                                 p_proj22(Ny)                
  !p: polarization   f: density		coup: magnetism coupling term
  complex*16                                   ::pt(2*Nm_o+1) = 0.0d0, &
                                                 ft(2*Nm_o+1) = 0.0d0, &
                                                 J_THZ_t = 0.0d0                 
  !pt: macroscopic polarization as a function of time ft: macroscopic density as a function of time
  complex*16                                   ::p_freq(N_freq) = 0.0d0, E_freq(N_freq) = 0.0d0,&
                                                 A_freq(N_freq) = 0.0d0, J_THZ_freq(N_freq) = 0.0d0, &
                                                 p1_freq(N_freq) = 0.0d0   
  !p_freq: Fourier transform of polarization   E_freq: Fourier transform of electrical field
  double precision                             ::freqgrid(N_freq) = 0.0d0, freqgrid1(N_freq) = 0.0d0, test(Ny) = 0.0d0, readin = 0.0d0, y_help(N_sec+1)
  !freqgrid: frequency grid used for Fourier transform   test: irrelevant with code
  character(80)                                ::list_file, list_file1, list_file2, list_file3, &
                                                 list_file4, list_file5
  !list_file: file name for output
  character(len=100)                           ::format_V, format_V1, format_V2, format_V3, format_V4     
  !format_V: format for output
  contains


  subroutine constant
  use N_of_grids
  !07/03/2017 creation
  !calculate y grid and step length for finder y grid
    integer                                      ::Ndo, Ndo_in   !number for do loop
!    dy = 0.08d0
!    ymax = dy*dble(Ny)
    ymax = 15.0d0
!    dy = ymax/dble(Ny) 
!    y_sec = Ny1*dy
y_sec = (/1.0d0,1.0d0,13.0d0/)
    do Ndo = 1, N_sec
      dy_sec(Ndo) = y_sec(Ndo)/Ny1(Ndo)
    end do

    Ny1_help(1) = 0
    y_help(1) = 0.0d0
    do Ndo = 1,N_sec
      Ny1_help((Ndo+1):(N_sec+1)) = Ny1_help((Ndo+1):(N_sec+1)) + Ny1(Ndo)
      y_help((Ndo+1):(N_sec+1)) = y_help((Ndo+1):(N_sec+1))+y_sec(Ndo)
    end do



    do Ndo = 1,N_sec
      do Ndo_in = (Ny1_help(Ndo)+1),Ny1_help(Ndo+1)
        y(Ndo_in) = y_help(Ndo)+dy_sec(Ndo)*(dble(Ndo_in - Ny1_help(Ndo)) - 0.50d0)           !y grid
      end do
    end do
    dy_fine = dy/dble(N_fine)             !y step for finer grid
    dy_fine_sec = dy_sec/dble(N_fine)             !y step for finer grid
do Ndo = 1,N_sec
    dy1((Ny1_help(Ndo)+1):Ny1_help(Ndo+1)) = dy_sec(Ndo)
!    write(*,*) (Ny1_help(Ndo)+1),Ny1_help(Ndo+1)
end do
!write(*,*) y

!    dipole = dipole/(1d0+10d0*y*y*E_B/Eg)
!----------------2 photon----------------------
!    A_freq_para = (Eg-0.4344975785796765d0*E_B)/hbar/2d0
!    A_freq_para = (Eg)/hbar/2d0
!----------------2 photon----------------------
!    A_freq_para = (Eg-0.07490971614016068d0*E_B)/hbar/4d0
!    A_freq_para = (Eg-0.1512538962722577d0*E_B)/hbar/3d0
!    A_excit = 1d4*(Eg/3d0/hbar)/A_freq_para

!    A_excit = 3.0d0/A_freq_para*1d8

!    A_excit = 1d1*(1d3)**(1/10d0*readin)*(Eg/3d0/hbar)/A_freq_para
    decay_m = 5.0d0*gamma;

!    open(unit = 100, file = 'fort.42', status = 'old', action = 'read')
!    read(100,*) readin
!----------------2 photon----------------------
!    A_freq_para = (Eg-(0.4d0 + 0.02d0*readin)*E_B)/hbar/2d0
!    A_freq_para = (Eg-(-1.9d0 - 0.1d0*readin)*E_B)/hbar/2d0
!    A_freq_para = (Eg-(5d0 - 15d0*readin/150d0)*E_B)/hbar/2d0
!----------------2 photon----------------------
!----------------3 photon----------------------
!   A_freq_para = (Eg-(0.4d0 - 0.02d0*readin)*E_B)/hbar/3d0
!    A_freq_para = (Eg-(-1.9d0 - 0.1d0*readin)*E_B)/hbar/3d0
!    A_freq_para = (Eg-(5d0 - 15d0*readin/150d0)*E_B)/hbar/3d0
!----------------3 photon----------------------
!    Nt = 5d3*scale1*(21d0-readin)
!    A_freq_para = (Eg-(3.0d0 - 0.1d0*readin)*E_B)/hbar/2d0
!   A_freq_para = (Eg-2.204334288969075d0*E_B)/hbar
   A_freq_para = Eg/hbar/2.0d0-15.0d0!+dble(readin)*1.5d0

!    A_freq_para = (Eg+40.0d0*E_B)/hbar/2d0
!    A_excit = 1.0d-2*(1.0d4)**(readin/40d0)*1d8/A_freq_para
!    decay_m =5.0d0*gamma*E_B*10.0d0**((readin - 1d0)/10.0d0*3d0)/E_B
!    gamma = (1.0d0+0.1d0*readin)/E_B
!    scale1 = 6.0d0*hbar/(gamma*E_B)
    scale1 = 2.8d-3/2.0d0/dsqrt(log(2.0d0))
    sigmat_A = 1.0d0/(gamma/2.0d0/pi*(E_B/hbar))
    tstart_A = 0.0d0
!    t_end = 6.0d0*hbar/(gamma*E_B)
    t_end = 6.0d0/(gamma/2.0d0/pi*(E_B/hbar))
    dt = pi/(4.0d3*E_B/hbar)
    Nt = int(t_end/dt)
!write(*,*) Nt
!    Nt = 10
    Nt_RWA = Nt

    delta_1s = A_freq_para - omega_1s
!    A_excit = 2d0*sqrt(sqrt(1d0/sigmat_A))/A_freq_para*1d8		!fix E*tau
    A_excit = 20.0d0/A_freq_para !mV/A
!    A_excit = 1.0d0*(13.5d0)**(1.0d0/30.0d0)*1d8/A_freq_para
!    dt = t_end/dble(Nt)  
!write(*,*) E_B/hbar

  end subroutine constant

end module constants
