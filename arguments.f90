  module constants  
  !07/03/2017 Creation
  !Initialization and calculation of all the constants and global variables 
  use N_of_grids
  implicit none

  double precision, parameter                  ::t_end = 1d0*scale1            
  !t_end: whole time
!  integer                                      ::m(2*Nm_o+1) = 0
!  !m: order of circular harmonics
  double precision, parameter                  ::ymax = 35.0d0, dy = ymax/dble(Ny) 
  !ymax: maximum value for y   dy: length of one step for y grid
  double precision, parameter                  ::sigmat = 0.3d0, tstart = -3d0,&
                                                 sigmat_A = 0.2d0, tstart_A = -1d0  
!                                                 sigmat_A = 0.2d0*scale1, tstart_A = -0.5d0*scale1  
  !sigmat: sigma for gaussian pulse(E(t) = exp(-t^2/sigmat^2))(ps)   tstart: starting point for
  !time (ps) 
  double precision, parameter                  ::pi = 4.0d0*datan(1.0d0)        
  !physical parameters
  double precision, parameter                  ::hbar = 0.6582119514d0, e = 1.60217662d-19  
  !physical parameters  unit of hbar: meV * ps
  double precision, parameter                  ::Ebind = 4.18d0, gamma = 1.56d0, Eg = 1.50d3
  !Ebind: binding energy(meV)   gamma: dephasing factor (meV) Eg: ground state energy
  double precision, parameter                  ::omega_1s = (0d0)/hbar, &    
                                                 A_freq_para = (Eg)/2d0/hbar
!                                                 A_freq_para = (4d0-4d0/9d0)*Ebind/hbar
!                                                 A_freq_para = omega_1s
!                                                dipole: unit: m*e(elementary charge)
  double precision, parameter                  ::dt = t_end/dble(Nt)           
  !length of one time step

!-------------------LAPACK-----------------------
  integer, parameter                           ::LDA = Ny, LDVL = Ny, LDVR = Ny
  integer, parameter                           ::LWMAX = 10000
  integer                                      ::INFO, LWORK
  double precision                             ::RWORK(2*Ny)
  complex*16                                   ::A( LDA, Ny ), VL( LDVL, Ny ), VR( LDVR, Ny ),&
                                                 W( Ny ), WORK( LWMAX ), VL1( LDVL, Ny ), &
                                                 VR1( LDVR, Ny ), W1( Ny ), A1( LDA, Ny ), &
                                                 VR_temp(LDVR, Ny), VL_temp(LDVL, Ny), &
                                                 VL2( LDVL, Ny ), VR2( LDVR, Ny ), &
                                                 W2( Ny ), A2( LDA, Ny )
  character                                    ::YES = 'V',NO = 'N'
!-------------------end LAPACK-----------------------
  double precision                             ::E_excit = 1.0d-3, shift = (Eg - omega_1s*hbar)/Ebind,A_excit = 1d6*(Eg/3d0/hbar)/A_freq_para
  !E_excit: excitation level(unit: binding energy)  shift: shift caused by rotation frame 
  !(unit: binding energy)
  !A_excit: unit: V*ps/m
  double precision                             ::coul_mat(Nm_o+1, Ny, Ny) = 0.0d0, &
                                                 Etemp(Ny), Etemp1(Ny), Emax, dipole(Ny) = 5d-10
  !coul_mat: Coulomb matrix(non-symmetric)(unit: binding energy) Et: electrical field 
  !for excitation
  double precision                             ::y(Ny)=0.0d0, delta_1s = A_freq_para - omega_1s
  !y: grid for y
  double precision                             ::y_fine(N_fine) = 0.0d0, dy_fine
  !y_fine: finer grid for removal of singularity   dy_fine: length of one step of finer grid for
  !removal of ringularity  f:density
  complex*16                                   ::p(2*Nm_o+1, Ny) = 0.0d0, ii = (0.0d0, 1.0d0), &
                                                 f(2*Nm_o+1, Ny) = 0.0d0, coup(Ny),&
                                                 p_proj(Nt_RWA, Ny), p_proj1(Nt_RWA, Ny), &
                                                 p_proj2(Nt_RWA, Ny), p_proj11(Nt_RWA, Ny), &
                                                 p_proj22(Nt_RWA, Ny)                
  !p: polarization   f: density		coup: magnetism coupling term
  complex*16                                   ::pt(2*Nm_o+1) = 0.0d0, &
                                                 ft(2*Nm_o+1) = 0.0d0, &
                                                 J_THZ_t = 0.0d0                 
  !pt: macroscopic polarization as a function of time ft: macroscopic density as a function of time
  complex*16                                   ::p_freq(N_freq) = 0.0d0, E_freq(N_freq) = 0.0d0,&
                                                 A_freq(N_freq) = 0.0d0, J_THZ_freq(N_freq) = 0.0d0, &
                                                 p1_freq(N_freq) = 0.0d0   
  !p_freq: Fourier transform of polarization   E_freq: Fourier transform of electrical field
  double precision                             ::freqgrid(N_freq) = 0.0d0, test(Nphi) = 0.0d0     
  !freqgrid: frequency grid used for Fourier transform   test: irrelevant with code
  character(80)                                ::list_file, list_file1, list_file2, list_file3, &
                                                 list_file4, list_file5
  !list_file: file name for output
  character(len=100)                           ::format_V, format_V1, format_V2, format_V3     
  !format_V: format for output
  contains


  subroutine constant
  !07/03/2017 creation
  !calculate y grid and step length for finder y grid
    integer                                      ::Ndo   !number for do loop
    do Ndo = 1,Ny
      y(Ndo) = dy*(dble(Ndo) - 0.50d0)           !y grid
    end do
    dy_fine = dy/dble(N_fine)             !y step for finer grid
    dipole = dipole/(1+y*y*Ebind/Eg)
  end subroutine constant

end module constants
