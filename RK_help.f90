  module RK_help
  !07/03/2017 Creation
  !All the subroutines and functions that Runge-Kutta needs, including calculation of E(t), time evolution of polarization and density, import frequency grids used for FFT
  use constants
  use N_of_grids
  implicit none
  contains


  subroutine readdata
  !07/03/2017 creation
  !import the grid of frequency grid Ben used, to make sure we are using the same frequency grid
  !08/20/2017 revision
  !generate a frequency grid from 1 to 9 binding energy
  integer, parameter                                        ::N_div = 1, N_sec = N_freq/N_div
  integer                                                   ::i1,i2
  !i1: do loop parameter
  !N_div: number of sections in frequency, e.g. omega, 2omega, 3omega... etc. 
  double precision                                          ::freq_bas(N_sec)
  !
  !--------set frequence grid, from 1 E_B to 9 E_B ---

  do i1 = 1,N_sec
    freq_bas(i1) =1d0/dble(N_freq)*(i1-1.0d0)-0.5d0
!    freq_bas(i1) =1d0/N_freq*i1-0.5d0/N_freq-0.5d0
!    freq_bas(i1) =1d0+0.01d0*i1-0.005d0
  end do
!-----------output the frequency grid---------------
!  freqgrid1 = 0
!  do i2 = 1,N_div
!    do i1 = 1,N_sec
!!      freqgrid(i1+(i2-1)*N_sec) = (A_freq_para*hbar*dble(i2)/4.0d0*3.0d0+A_freq_para*hbar/2.0d0*freq_bas(i1))/hbar
!      freqgrid(i1+(i2-1)*N_sec) = (A_freq_para*hbar*dble(i2)+freq_bas(i1)*10d0*E_B)/hbar
!    end do
!  end do
  freqgrid = (Eg+freq_bas*20.0d0*E_B)/hbar
!  freqgrid = (1.5d0*Eg+freqgrid*3d0*Eg)/hbar
  write(list_file, '(A)') 'freqgrid.dat'           !p(t)
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) freqgrid(i1)*hbar/Eg
  end do
  close(700)
  write(list_file, '(A)') 'freqgrid1.dat'           !p(t)
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) freqgrid1(i1)*hbar/Eg
  end do
  close(700)
!------------end output-----------------------------
  end subroutine readdata

  complex*16 function Etime(tvia)
  !07/03/2017 creation
  !function for calculating electrical field(unit: binding energy)
    double precision                             ::tvia
    double precision                             ::t1
    t1 = tvia*dt+tstart_A -3.0d0*sigmat_A
    !tvia: time step(converted into double precision)
!    if (t1 - tstart_A < 80.0d0*sigmat_A) then
!        Etime = A_excit*1.0d3*(sin(A_freq_para*t1))*A_freq_para!meV non RWA, 1d3 for eV to meV        
!    else
!	Etime = 0.0d0
!    end if  
    Etime = A_excit*exp(-t1*t1/(sigmat_A*sigmat_A)) * (A_freq_para*(exp(ii*delta_1s*t1)*exp(ii*2.0d0*omega_1s*t1)-exp(-ii*delta_1s*t1))/(2.0d0*ii) + 2.0d0*t1/(sigmat_A*sigmat_A) * (exp(-ii*delta_1s*t1) + exp(ii*delta_1s*t1)*exp(ii*2.0d0*omega_1s*t1))/2.0d0)!meV non RWA, 1d3 for eV to meV 
!    Etime = A_excit*1.0d3*(sin(A_freq_para*t1))!meV non RWA, 1d3 for eV to meV 
!    Etime = A_excit*E_B*exp(-t1*t1/(sigmat_A*sigmat_A)) * (A_freq_para*sin(A_freq_para*t1))/exp(ii*omega_1s*t1)!no E_B non RWA just sinusoidal part
!    Etime = E_excit*exp(-t1*t1/(sigmat*sigmat)) * (1d0 + exp(-ii*2.0d0*omega_1s*t1))!no E_B   RWA
!    Etime = A_excit*exp(-(t1)*(t1)/(sigmat_A*sigmat_A))*dcos(A_freq_para*t1) !unit: mV/A
  end function Etime

  complex*16 function Atime(tvia)
  !09/03/2017 creation
  !function for calculating magnetic field(unit: binding energy)
    double precision                             ::tvia  
    double precision                             ::t1
    t1 = tvia*dt+tstart_A -3.0d0*sigmat_A                  
    !tvia: time step(converted into double precision)
    Atime = A_excit*exp(-t1*t1/(sigmat_A*sigmat_A))*(cos(A_freq_para*t1))        ! unit:mV*ps/A           
!    Atime = A_excit*exp(-(t1)*(t1)/(sigmat_A*sigmat_A))*cos(A_freq_para*t1/hbar)*E_B*2.0d0
!    Atime = A_excit*exp(-(t1)*(t1)/(sigmat_A*sigmat_A))*(exp(ii*A_freq_para*t1/hbar)+exp(-ii*A_freq_para*t1/hbar))*E_B        ! unit:meV  
!    if (t1 - tstart_A < 80.0d0*sigmat_A) then
!        Atime = A_excit*(cos(A_freq_para*t1))        ! unit:V*ps/m         
!    else
!	Atime = 0.0d0
!    end if     
  end function Atime


  subroutine RHS(nt_via, f_via, p_via, decay_via, &
                 p_out, f_out, decay_out, kE_out, kPfreq_out, kP1freq_out, kA_out, kJ_out)
  !08/25/2017 creation
  !calculate the right hand side of RK
    double precision, intent(in)                           ::nt_via   
    !time step(converted into double precision)
    complex*16, intent(in)                                 ::p_via(2*Nm_o+1, Ny)  
    !polarization of last time step
    complex*16, intent(in)                                 ::f_via(2*Nm_o+1, Ny)   
    !density of last time step
    complex*16, intent(in)                                 ::decay_via(2*Nm_o+1, Ny)   
    complex*16, intent(out)                                ::p_out(2*Nm_o+1, Ny)  
    !polarization output
    complex*16, intent(out)                                ::f_out(2*Nm_o+1, Ny) 
    !density output
    complex*16, intent(out)                                ::decay_out(2*Nm_o+1, Ny) 
    complex*16, intent(out)                                ::kE_out(N_freq)
    complex*16, intent(out)                                ::kPfreq_out(N_freq)
    complex*16, intent(out)                                ::kP1freq_out(N_freq)
    complex*16, intent(out)                                ::kA_out(N_freq)
    complex*16, intent(out)                                ::kJ_out(N_freq)
    !kE_out, kPfreq_out, kA_out, kJ_out: output parameters for Runge-Kutta methods
    integer                                                ::Ndo_m, Ndo_m1, Ndo_phi
    !Ndo_m, Ndo_m1: arguments for do loops, indices of m and m'
    complex*16                                             ::p_sum_part(Ny), &
                                                             f_sum_part(Ny), &
                                                             pp_sum(Ny), &
                                                             pp_sum_plus(Ny), &
                                                             pf_sum(Ny), &
                                                             fp_sum(Ny), &
                                                             p_sum_part_m(Ny), &
                                                             decayf_sum(Ny), &
                                                             fdecay_sum(Ny), &
                                                             decay_sum_part(Ny), &
                                                             decay_sum_part_m(Ny)
    double precision                                       ::t1
    !help parameters to compute all the summations
    !p_sum_part: \Sigma V_{k'-k}P_{k'}
    !f_sum_part: \Sigma V_{k'-k}f_{k'}
    !pp_sum: \Sigma_{m'} P^*_{k, -m+m'} \Sigma_{k'}V^{m'}_{k,k'} P_{k', m'}
    !pp_sum_plus: \Sigma_{m'} P_{k, m+m'} \Sigma_{k'}V^{-m'}_{k,k'} P^*_{k', m'}
    !pf_sum: \Sigma_{m'} P^*_{k, m-m'} \Sigma_{k'}V^{m'}_{k,k'} f_{k', m'}
    !fp_sum: \Sigma_{m'} f^*_{k, m-m'} \Sigma_{k'}V^{m'}_{k,k'} P_{k', m'}
    !p_sum_part_m: \Sigma V^m_{k, k'}P^m_{k'}
    complex*16                                             ::pt1, J_THZ_t1
    !pt1:\Sigma_k P_{k}
    !J_THZ_t1: \Sigma_k J_{k} f_k

    pt1 = 0.0d0
    !Initialization of P(t)
    do Ndo_m = 1, 2*Nm_o+1
      p_sum_part = 0.0d0
      pp_sum = 0.0d0
      pp_sum_plus = 0.0d0
      pf_sum = 0.0d0
      fp_sum = 0.0d0
      p_sum_part_m = 0.0d0
      f_sum_part = 0.0d0
      decayf_sum = 0.0d0
      fdecay_sum = 0.0d0
      decay_sum_part = 0.0d0
      decay_sum_part_m = 0.0d0
      do Ndo_m1 = 1, 2*Nm_o+1

        f_sum_part = matmul(f_via(Ndo_m1, :), coul_mat_f(abs(Ndo_m1-Nm_o-1)+1, :, :))                                                                               
        p_sum_part = matmul(p_via(Ndo_m1, :), coul_mat_ch(abs(Ndo_m1-Nm_o-1)+1, :, :))
        decay_sum_part = matmul(decay_via(Ndo_m1, :), coul_mat_ch(abs(Ndo_m1-Nm_o-1)+1, :, :))                        
        if(abs(Ndo_m-Ndo_m1)<Nm_o+1) then
          fp_sum = fp_sum+f_via(Ndo_m-Ndo_m1+(Nm_o+1), :)*p_sum_part
          pf_sum = pf_sum+p_via(Ndo_m-Ndo_m1+(Nm_o+1), :)*f_sum_part
          fdecay_sum = fdecay_sum+f_via(Ndo_m-Ndo_m1+(Nm_o+1), :)*decay_sum_part
          decayf_sum = decayf_sum+decay_via(Ndo_m-Ndo_m1+(Nm_o+1), :)*f_sum_part
          pp_sum = pp_sum+conjg(p_via(-Ndo_m+Ndo_m1+(Nm_o+1), :))*p_sum_part
        end if
        if(abs(Ndo_m+Ndo_m1 - 2*(Nm_o+1))<Nm_o+1) then
          pp_sum_plus = pp_sum_plus+p_via(Ndo_m+Ndo_m1-(Nm_o+1), :)*conjg(p_sum_part)
        end if
      end do

      !magnetism coupling, P^{m-1}+P^{m+1}, and truncate at Nm_o, -Nm_o
      if(abs(Ndo_m-Nm_o-1)<Nm_o) then
        coup = (p_via(Ndo_m-1, :)+p_via(Ndo_m+1, :))/2.0d0
        else if(Ndo_m == 1) then
        coup = (p_via(Ndo_m+1, :))/2.0d0
        else if(Ndo_m == 2*Nm_o+1) then
        coup = (p_via(Ndo_m-1, :))/2.0d0
      end if
!      coup = 0d0
      p_sum_part_m = matmul(p_via(Ndo_m, :), coul_mat_ch(abs(Ndo_m-Nm_o-1)+1, :, :))
      p_out(Ndo_m, :) = -ii*(y*y*p_via(Ndo_m, :) - 2.0d0*pf_sum - 2.0d0*a_B/hbar*y*Atime(nt_via)*coup&
                        +shift*p_via(Ndo_m, :)&
                        -((abs(dble(Ndo_m == (Nm_o+1)))-2.0d0*f_via(Ndo_m, :)) *Etime(nt_via)*dipole/E_B+&
                       (p_sum_part_m - 2.0d0*fp_sum) )- ii*decay_m*decay_via(Ndo_m, :) )/hbar*dt*E_B
!      p_out(Ndo_m, :) = -(0.0d0,1.0d0)*(y*y*p_via(Ndo_m, :) - 2.0d0*pf_sum - 2.0d0*a_B/hbar*1.0d3*y*Atime(nt_via)*coup+&
!                        shift*p_via(Ndo_m, :) - (0.0d0,1.0d0) * gamma * p_via(Ndo_m, :) -&
!                        ((abs(dble(Ndo_m == (Nm_o+1)))-2.0d0* f_via(Ndo_m, :)) *Etime(nt_via)*dipole/E_B+&
!                        (p_sum_part_m - 2.0d0*fp_sum) ) )/hbar*dt*E_B
!      f_out = 0.0d0
      f_out(Ndo_m, :) = -ii*(conjg(Etime(nt_via)*dipole/E_B)*p_via(Ndo_m, :) - Etime(nt_via)*dipole/E_B&
                        *conjg(p_via((2*Nm_o+2)-Ndo_m, :)) &
                        +(pp_sum_plus-pp_sum)- ii *(dble(Ndo_m==(Nm_o+1))+1.0d0)*gamma &
                        * f_via(Ndo_m, :))/hbar*dt*E_B
      decay_sum_part_m = matmul(decay_via(Ndo_m, :), coul_mat_ch(abs(Ndo_m-Nm_o-1)+1, :, :))
      decay_out(Ndo_m, :) = -ii*(y*y*decay_via(Ndo_m, :) - ii*decay_m* decay_via(Ndo_m,:) - 2.0d0*decayf_sum + &
                        shift*decay_via(Ndo_m, :) + ii * gamma * p_via(Ndo_m, :) -&
                        ((decay_sum_part_m-2.0d0*fdecay_sum) ) )/hbar*dt*E_B
    end do
!      write(807, format_V) -real((conjg(Etime(nt_via)*dipole/E_B)*p_via(Nm_o+1, :) - Etime(nt_via)*dipole/E_B&
!                        *conjg(p_via(Nm_o+1, :)))&
!                        /hbar*dt*ii*E_B)
!      write(808, format_V) real(f_out(Nm_o+1,:))
    t1 = dt*nt_via+tstart_A
    J_THZ_t1 = 0d0
!    J_THZ_t1 = dy1*sum(y*y*(f_via(Nm_o, :)+f_via(Nm_o+2, :)))/(2.0d0*pi)
    kE_out = dt * Etime(nt_via) * exp(ii*t1*freqgrid)  
    kA_out = dt * Atime(nt_via) * exp(ii*t1*freqgrid) 
    pt1 = sum(dy1*y*p_via(Nm_o+1, :)*dipole)/(2.0d0*pi)/a_B/a_B   !calculate macroscopic polarization unit: e/m
    kPfreq_out = dt * pt1 * exp(ii*freqgrid*t1)
    kJ_out = dt * J_THZ_t1 * exp(ii*freqgrid*t1)
    kP1freq_out = dt * p_proj2(1)/(2.0d0*pi) * exp(ii*freqgrid*t1)
  end subroutine RHS
  end module RK_help
