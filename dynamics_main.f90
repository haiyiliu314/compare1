program main
  !07/03/2017 creation
  !main program to calculate dynamics of SBE
  !result: susceptibility
  !zero-density
  !07/06/2017 Revision   
  !remove RWA
  !07/07/2017 Upgrade   
  !get rid of arrays with dimension(1,Nt), use elements to increase the speed
  !08/23/2017 Revision
  !calculate the circular harmonics expansion
  use constants
  use RK_help
  use N_of_grids
  use coul_mat_module
  use RK_module
  use params
  implicit none
  integer                                                       ::Ndo, Ndo_m, i2, Ndo_in, i3, i4
  complex*16                                                    ::osc_str(Ny)
  call constant  
  call coul_matrix
  call readdata

  write(format_V1, '(A12, I8, A18)')   '(SE48.32e3, ', 1, '(", ",SE48.32e3))'
  write(format_V2, '(A12, I8, A18)')   '(SE48.32e3, ', 2*Nm_o+1, '(", ",SE48.32e3))'
  write(format_V3, '(A12, I8, A18)')   '(SE48.32e3, ', 2, '(", ",SE48.32e3))'
  write(format_V4, '(A12, I8, A18)')   '(SE48.32e3, ', Ny, '(", ",SE48.32e3))'


!  write(list_file, '(A)') 'f_end.dat'             !f
!  open(unit=709,file=list_file)
!  do i4 = 10:10
!  A_excit = (1+dble(i4)*0.1)*1d5*(Eg/3.0d0/hbar)/A_freq_para





  write(list_file, '(A)') 'pk2.dat'           !pk
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE48.32e3, ', Ny, '(", ",SE48.32e3))'
  write(list_file, '(A)') 'p_r.dat'             !p
  open(unit=701,file=list_file)
  write(list_file2, '(A)') 'Et.dat'               !E(t)
  open(unit=702,file=list_file2)
  write(list_file, '(A)') 'pt.dat'           !p(t)
  open(unit=703,file=list_file)
  write(list_file, '(A)') 'p1t.dat'           !p(t)
  open(unit=711,file=list_file)
  write(list_file, '(A)') 'p_i.dat'           !p(t)
  open(unit=710,file=list_file)
  write(list_file, '(A)') 'ft.dat'           !f(t)
  open(unit=704,file=list_file)
  write(list_file, '(A)') 'At.dat'           !A(t)
  open(unit=705,file=list_file)
  write(list_file, '(A)') 'Jt.dat'           !J(t)
  open(unit=706,file=list_file)
  write(list_file, '(A)') 'f_r.dat'             !f
  open(unit=707,file=list_file)
  write(list_file, '(A)') 'f_i.dat'             !f
  open(unit=708,file=list_file)
  write(list_file, '(A)') 'p_proj.dat'       !eigen values
  open(unit=800,file=list_file)
  write(list_file, '(A)') 'p_proj1.dat'       !eigen values
  open(unit=801,file=list_file)
  write(list_file, '(A)') 'p_proj2.dat'       !eigen values
  open(unit=802,file=list_file)
  write(list_file, '(A)') 'decay_r.dat'             !f
  open(unit=807,file=list_file)
  write(list_file, '(A)') 'decay_i.dat'             !f
  open(unit=808,file=list_file)
  i2 = 1
  i3 = 1


  do Ndo = 1, Nt
    
    !calculate p(0)
    do Ndo_in = 1, 2*Nm_o+1
      ft(Ndo_in) = dy*sum(y*f(Ndo_in, :))/(2.0d0*pi)   !calculate macroscopic density
      pt(Ndo_in) = dy*sum(y*p(Ndo_in, :))/(2.0d0*pi)/(0.125d0)**2   !calculate macroscopic polarization
    end do
    J_THZ_t = dy*sum(y*y*(f(Nm_o, :)+ f(Nm_o+2, :)))/(4.0d0*pi)   !calculate macroscopic density
    if(Ndo == i2*Nt/Nt_RWA) then
      write(703, format_V2) abs(pt)
      write(704, format_V2) abs(ft)
      write(711, format_V3) abs(dy*sum(y*y*(p(Nm_o, :)+ p(Nm_o+2, :)))/(4.0d0*pi))
    end if
    call RK(E_freq, p_freq, p1_freq, f, p, decay, real(Ndo,8), A_freq, J_THZ_freq, &              !input
            E_freq, p_freq, p1_freq, f, p, decay, A_freq, J_THZ_freq)                           !output
    if(Ndo == i2*Nt/Nt_RWA) then
     
      write(702, format_V1) real(Etime(real(Ndo,8)))
      write(705, format_V3) real(Atime(real(Ndo,8))), aimag(Atime(real(Ndo,8)))
      write(706, format_V3) real(J_THZ_t), aimag(J_THZ_t)
      i2 = i2+1
    end if
    if(Ndo == i3*Nt/Nt_RWA) then
      do i1 = 1, 2*Nm_o+1
        write(701, format_V) real(p(i1, :))
        write(710, format_V) aimag(p(i1, :))
        write(707, format_V) real(f(i1, :))
        write(708, format_V) aimag(f(i1, :))
        write(807, format_V) real(decay(i1, :))
        write(808, format_V) aimag(decay(i1, :))
      end do
      i3 = i3+1
    end if
  end do
  close(700)
  close(701)
  close(704)
  close(703)
  close(705)
  close(706)
  close(707)
  close(708)
  close(711)
  close(800)
  close(801)
  close(802)
  close(807)
  close(808)
!  end do
!-------------------output--------------------
  write(list_file, '(A)') 'coulomb.dat'       !coul_mat
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE48.32e3, ', Ny, '(", ",SE48.32e3))'
  do i1 = 1, Nm_o+1
    do i2 = 1, Ny
      write(700, format_V) real(coul_mat(i1, i2, :))
    end do
  end do
  close(700)
  write(list_file, '(A)') 'p_freq.dat'           !p_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE48.32e3, ', N_freq, '(", ",SE48.32e3))'
  do i1 = 1, N_freq
    write(700, format_V) p_freq(i1)
  end do
  close(700)

  write(list_file, '(A)') 'p1_freq.dat'           !p_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE48.32e3, ', N_freq, '(", ",SE48.32e3))'
  do i1 = 1, N_freq
    write(700, format_V) p1_freq(i1)
  end do
  close(700)
    
  write(list_file, '(A)') 'E_freq.dat'           !E_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE48.32e3, ', N_freq, '(", ",SE48.32e3))'
  do i1 = 1, N_freq
    write(700, format_V) E_freq(i1)
  end do
  close(700)

  write(list_file, '(A)') 'A_freq.dat'           !E_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE48.32e3, ', N_freq, '(", ",SE48.32e3))'
  do i1 = 1, N_freq
    write(700, format_V) A_freq(i1)
  end do
  close(700)

  write(list_file, '(A)') 'J_freq.dat'           !E_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE48.32e3, ', N_freq, '(", ",SE48.32e3))'
  do i1 = 1, N_freq
    write(700, format_V) J_THZ_freq(i1)
  end do
  close(700)

  write(list_file, '(A)') 'y.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE48.32e3, ', Ny, '(", ",SE48.32e3))'
      write(700, format_V) y
  close(700)

  write(list_file, '(A)') 'dipole.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE48.32e3, ', Ny, '(", ",SE48.32e3))'
      write(700, format_V) dipole
  close(700)


  write(list_file, '(A)') 'para_mat.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE48.32e3, ', 1, '(", ",SE48.32e3))'  
    write(700, format_V) t_end		!1
    write(700, format_V) dt*1d6   		!2
    write(700, *) Nt		!3
    write(700, format_V) ymax		!4
    write(700, *) Ny		!5
    write(700, format_V) A_excit*A_freq_para*1d-8		!6
    write(700, *) gamma		!7
    write(700, *) 2.0d0*pi/A_freq_para*1000.0d0		!8
    write(700, *) Nt_RWA		!9
    write(700, *) tstart_A		!10
    write(700, *) sigmat_A		!11
    write(700, *) Nm_o		!12
    write(700, *) N_freq		!13
    write(700, *) delta_1s		!14
    write(700, *) decay_m		!15
    write(700, *) gamma		!16
  close(700)

!----------------end output------------------
  call output_params
end program main
