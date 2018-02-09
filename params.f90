module params
!07/03/2017 creation
!output main parameters that used in the run
  use constants
  use RK_help
  use N_of_grids
  implicit none
  contains
  subroutine output_params
    write(list_file1, '(A)') 'params.dat'    
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', 1, '(", ",SE24.16e3))'   
    open(unit=600,file=list_file1)  
    write(600, *) 'Time length(ps)' 
    write(600, format_V) t_end
    write(600, *) 'dt(as)' 
    write(600, format_V) dt*1d6   
    write(600, *) 'Number of time steps' 
    write(600, *) Nt
    write(600, *) 'maximum value of y' 
    write(600, format_V) ymax
    write(600, *) 'Peak of electrical field(MV/cm)' 
    write(600, format_V) A_excit*Ebind*1.129572800243224d3/50d0
    write(600, *) 'dephasing(meV)'
    write(600, *) gamma
    write(600, *) 'time for one cycle(fs)'
    write(600, *) 2d0*pi/A_freq_para
  end subroutine output_params
end module params
