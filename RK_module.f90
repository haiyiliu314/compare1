module RK_module
!07/03/2017 creation
!Runge-Kutta for: FFT of electrical field & polarization; time evolution of polarization and density
  use constants
  use RK_help
  use N_of_grids
  implicit none
  contains
  subroutine RK(E_in, Pfreq_in, P1freq_in, f_in, p_in, decay_in, n, A_in, J_THZ_in, &
                E_out, Pfreq_out, P1freq_out, f_out, p_out, decay_out, A_out, J_THZ_out)
  !07/03/2017 creation
  !Runga_Kutta subroutine
    double precision                                    ::n                   !time step
    complex*16, intent(in)                              ::f_in(2*Nm_o+1, Ny)  !input density
    complex*16, intent(in)                              ::p_in(2*Nm_o+1, Ny)  !input polarization
    complex*16, intent(in)                              ::decay_in(2*Nm_o+1, Ny)
    complex*16, intent(in)                              ::E_in(N_freq)        !input FT of electrical field
    complex*16, intent(in)                              ::A_in(N_freq)        !input FT of vectorial potential
    complex*16, intent(in)                              ::J_THZ_in(N_freq)    !input FT of THz current
    complex*16, intent(in)                              ::Pfreq_in(N_freq)    !input FT of polarization
    complex*16, intent(in)                              ::P1freq_in(N_freq)    !input FT of polarization
    complex*16, intent(out)                             ::f_out(2*Nm_o+1, Ny) !output density
    complex*16, intent(out)                             ::p_out(2*Nm_o+1, Ny) !output polarization
    complex*16, intent(out)                             ::decay_out(2*Nm_o+1, Ny)
    complex*16, intent(out)                             ::E_out(N_freq)       !output FT of electrical field
    complex*16, intent(out)                             ::A_out(N_freq)       !output FT of vectorial potential
    complex*16, intent(out)                             ::J_THZ_out(N_freq)   !output FT of THz current
    complex*16, intent(out)                             ::Pfreq_out(N_freq)   !output FT of polarization
    complex*16, intent(out)                             ::P1freq_out(N_freq)   !output FT of P_p
    complex*16, dimension(2*Nm_o+1, Ny)                 ::kf1, kf2, kf3, kf4, kf5, kf6  
    complex*16, dimension(2*Nm_o+1, Ny)                 ::kp1, kp2, kp3, kp4, kp5, kp6
    complex*16, dimension(2*Nm_o+1, Ny)                 ::kdecay1, kdecay2, kdecay3, kdecay4, kdecay5, kdecay6
    complex*16, dimension(N_freq)                       ::kE1, kE2, kE3, kE4, kE5, kE6, &
                                                          kA1, kA2, kA3, kA4, kA5, kA6, &
                                                          kJ1, kJ2, kJ3, kJ4, kJ5, kJ6
    complex*16, dimension(N_freq)                       ::kPfreq1, kPfreq2, kPfreq3, kPfreq4,kPfreq5, kPfreq6,& 
                                                          kP1freq1, kP1freq2, kP1freq3, kP1freq4, kP1freq5, kP1freq6

    !k(...)1: the first step of Runge-Kutta(RK), and k(...)2 is the second step, etc. 
    n = n-1d0         !to fix the shift of one time step
    
    !calculate k1
    call RHS(n, f_in, p_in, decay_in, &
             kp1, kf1, kdecay1, kE1, kPfreq1, kP1freq1, kA1, kJ1)

    !calculate k2
    call RHS(n+0.250d0, f_in + kf1/4.0d0, p_in + kp1/4.0d0, decay_in + kdecay1/4.0d0, &
             kp2, kf2, kdecay2, kE2, kPfreq2, kP1freq2, kA2, kJ2)

    !calculate k3
    call RHS(n+3d0/8d0, f_in + 3d0/32.0d0*kf1+9d0/32.0d0*kf2, p_in + 3d0/32.0d0*kp1+9d0/32.0d0*kp2, decay_in + 3d0/32.0d0*kdecay1+9d0/32.0d0*kdecay2, &
             kp3, kf3, kdecay3, kE3, kPfreq3, kP1freq3, kA3, kJ3)

    !calculate k4
    call RHS(n+12d0/13d0, f_in + 1932d0/2197d0*kf1-7200d0/2197d0*kf2+7296d0/2197d0*kf3, p_in + 1932d0/2197d0*kp1-7200d0/2197d0*kp2+7296d0/2197d0*kp3, decay_in + 1932d0/2197d0*kdecay1-7200d0/2197d0*kdecay2+7296d0/2197d0*kdecay3, &
             kp4, kf4, kdecay4, kE4, kPfreq4, kP1freq4, kA4, kJ4)
    call RHS(n+1d0, f_in+439d0/216d0*kf1-8d0*kf2+3680d0/513d0*kf3-845d0/4104d0*kf4, p_in+439d0/216d0*kp1-8d0*kp2+3680d0/513d0*kp3-845d0/4104d0*kp4, decay_in+439d0/216d0*kdecay1-8d0*kdecay2+3680d0/513d0*kdecay3-845d0/4104d0*kdecay4, &
             kp5, kf5, kdecay5, kE5, kPfreq5, kP1freq5, kA5, kJ5)    
    call RHS(n+1d0/2d0, f_in-8d0/27d0*kf1+2d0*kf2-3544d0/2565d0*kf3+1859d0/4104d0*kf4-11d0/40d0*kf5, p_in-8d0/27d0*kp1+2d0*kp2-3544d0/2565d0*kp3+1859d0/4104d0*kp4-11d0/40d0*kp5, decay_in-8d0/27d0*kdecay1+2d0*kdecay2-3544d0/2565d0*kdecay3+1859d0/4104d0*kdecay4-11d0/40d0*kdecay5, &
             kp6, kf6, kdecay6, kE6, kPfreq6, kP1freq6, kA6, kJ6)
    !final step for Runge-Kutta
    p_out = p_in + kp1*16d0/135d0 +kp3*6656d0/12825d0 + kp4*28561d0/56430d0 + kp5*(-9d0/50d0) + kp6*2d0/55d0
    decay_out = decay_in + kdecay1*16d0/135d0 +kdecay3*6656d0/12825d0 + kdecay4*28561d0/56430d0 + kdecay5*(-9d0/50d0) + kdecay6*2d0/55d0
    f_out = f_in + kf1*16d0/135d0 +kf3*6656d0/12825d0 + kf4*28561d0/56430d0 + kf5*(-9d0/50d0)  + kf6*2d0/55d0
    E_out = E_in + kE1*16d0/135d0 +kE3*6656d0/12825d0 + kE4*28561d0/56430d0 + kE5*(-9d0/50d0) + kE6*2d0/55d0
    A_out = A_in + kA1*16d0/135d0 +kA3*6656d0/12825d0 + kA4*28561d0/56430d0 + kA5*(-9d0/50d0) + kA6*2d0/55d0
    J_THZ_out = J_THZ_in + kJ1*16d0/135d0 +kJ3*6656d0/12825d0 + kJ4*28561d0/56430d0 + kJ5*(-9d0/50d0) + kJ6*2d0/55d0
    Pfreq_out = Pfreq_in + kPfreq1*16d0/135d0 +kPfreq3*6656d0/12825d0 + kPfreq4*28561d0/56430d0 + kPfreq5*(-9d0/50d0) + kPfreq6*2d0/55d0
    P1freq_out = P1freq_in + kP1freq1*16d0/135d0 +kP1freq3*6656d0/12825d0 + kP1freq4*28561d0/56430d0 + kP1freq5*(-9d0/50d0) + kP1freq6*2d0/55d0
    end subroutine RK
end module RK_module
