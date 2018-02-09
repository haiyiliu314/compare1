module RK_module
!07/03/2017 creation
!Runge-Kutta for: FFT of electrical field & polarization; time evolution of polarization and density
  use constants
  use RK_help
  use N_of_grids
  implicit none
  contains
  subroutine RK(E_in, Pfreq_in, f_in, p_in, n, A_in, J_THZ_in, &
                E_out, Pfreq_out, f_out, p_out, A_out, J_THZ_out)
  !07/03/2017 creation
  !Runga_Kutta subroutine
    double precision                                    ::n                   !time step
    complex*16, intent(in)                              ::f_in(2*Nm_o+1, Ny)  !input density
    complex*16, intent(in)                              ::p_in(2*Nm_o+1, Ny)  !input polarization
    complex*16, intent(in)                              ::E_in(N_freq)        !input FT of electrical field
    complex*16, intent(in)                              ::A_in(N_freq)        !input FT of vectorial potential
    complex*16, intent(in)                              ::J_THZ_in(N_freq)    !input FT of THz current
    complex*16, intent(in)                              ::Pfreq_in(N_freq)    !input FT of polarization
    complex*16, intent(out)                             ::f_out(2*Nm_o+1, Ny) !output density
    complex*16, intent(out)                             ::p_out(2*Nm_o+1, Ny) !output polarization
    complex*16, intent(out)                             ::E_out(N_freq)       !output FT of electrical field
    complex*16, intent(out)                             ::A_out(N_freq)       !output FT of vectorial potential
    complex*16, intent(out)                             ::J_THZ_out(N_freq)   !output FT of THz current
    complex*16, intent(out)                             ::Pfreq_out(N_freq)   !output FT of polarization
    complex*16, dimension(2*Nm_o+1, Ny)                 ::kf1, kf2, kf3, kf4  
    complex*16, dimension(2*Nm_o+1, Ny)                 ::kp1, kp2, kp3, kp4
    complex*16, dimension(N_freq)                       ::kE1, kE2, kE3, kE4, &
                                                          kA1, kA2, kA3, kA4, &
                                                          kJ1, kJ2, kJ3, kJ4 
    complex*16, dimension(N_freq)                       ::kPfreq1, kPfreq2, kPfreq3, kPfreq4,& 
                                                          kPfreq1_hlp, kPfreq2_hlp, kPfreq3_hlp, &
                                                          kPfreq4_hlp 
    !k(...)1: the first step of Runge-Kutta(RK), and k(...)2 is the second step, etc. 
    n = n-1d0         !to fix the shift of one time step
    
    !calculate k1
    call RHS(n, f_in, p_in, &
             kp1, kf1, kE1, kPfreq1, kA1, kJ1)

    !calculate k2
    call RHS(n+0.50d0, f_in + kf1/2.0d0, p_in + kp1/2.0d0, &
             kp2, kf2, kE2, kPfreq2, kA2, kJ2)

    !calculate k3
    call RHS(n+0.50d0, f_in + kf2/2.0d0, p_in + kp2/2.0d0, &
             kp3, kf3, kE3, kPfreq3, kA3, kJ3)

    !calculate k4
    call RHS(n+1.0d0, f_in + kf3, p_in + kp3, &
             kp4, kf4, kE4, kPfreq4, kA4, kJ4)

    !final step for Runge-Kutta
    p_out = p_in + kp1/6.0d0 +kp2/3.0d0 + kp3/3.0d0 + kp4/6.0d0
    f_out = f_in + kf1/6.0d0 +kf2/3.0d0 + kf3/3.0d0 + kf4/6.0d0 
    E_out = E_in + kE1/6.0d0 +kE2/3.0d0 + kE3/3.0d0 + kE4/6.0d0
    A_out = A_in + kA1/6.0d0 +kA2/3.0d0 + kA3/3.0d0 + kA4/6.0d0
    J_THZ_out = J_THZ_in + kJ1/6.0d0 +kJ2/3.0d0 + kJ3/3.0d0 + kJ4/6.0d0
    Pfreq_out = Pfreq_in + kPfreq1/6.0d0 +kPfreq2/3.0d0 + kPfreq3/3.0d0 + kPfreq4/6.0d0

    end subroutine RK
end module RK_module
