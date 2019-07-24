module params
!07/03/2017 creation
!output main parameters that used in the run
  use constants
  use RK_help
  use N_of_grids
  implicit none
  contains
  subroutine output_params
    write(list_file1, '(A)') 'params.txt'    
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', 1, '(", ",SE24.16e3))'   
    open(unit=600,file=list_file1)  
    write(600, *) '============Constants================'
    write(600, *) 'pi' 
    write(600, format_V) pi
    write(600, *) 'LWORK' 
    write(600, format_V) LWORK1
    write(600, *) 'hbar(meV*ps)' 
    write(600, format_V) hbar
    write(600, *) 'elementary charge' 
    write(600, format_V) e
    write(600, *) 'fine-structure constant' 
    write(600, format_V) f_s
    write(600, *) 'speed of light in vacuum(A/ps)' 
    write(600, format_V) c_vac*1.0d10
    write(600, *) '============Material parameters================'
    write(600, *) 'E_B(meV)' 
    write(600, format_V) E_B
    write(600, *) 'E_gap(meV)' 
    write(600, format_V) Eg
    write(600, *) 'dipole(e*A)' 
    write(600, format_V) dipole(1)
    write(600, *) 'a_B(A)' 
    write(600, format_V) a_B
    write(600, *) 'dephasing(meV)' 
    write(600, format_V) gamma*E_B
    write(600, *) 'relative permittivity' 
    write(600, format_V) eps_r
    write(600, *) '============Pulse properties================'
    write(600, *) 'E_0(kV/cm)' 
    write(600, format_V) A_excit*1.0d2*A_freq_para
    write(600, *) 'pulse duration tau(ps)' 
    write(600, format_V) sigmat_A
    write(600, *) 'time for one cycle(ps)'
    write(600, *) 2d0*pi/A_freq_para
    write(600, *) '============Parameters for calculation================'
    write(600, *) 'starting time(ps)' 
    write(600, format_V) tstart_A
    write(600, *) 'Time length(ps)' 
    write(600, format_V) t_end
    write(600, *) 'tau(ps)' 
    write(600, format_V) sigmat_A
    write(600, *) 'dt(ps)' 
    write(600, format_V) dt  
    write(600, *) 'Number of time steps' 
    write(600, *) Nt
    write(600, *) 'maximum value of y' 
    write(600, format_V) ymax
    write(600, *) 'number of y grid' 
    write(600, *) Ny
    write(600, *) 'number of y grid in removing singularity' 
    write(600, *) N_fine


  end subroutine output_params
end module params
