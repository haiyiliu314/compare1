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
  double precision       :: AA(Ny,Ny)
  call constant
  call coul_matrix
  call readdata

  write(format_V1, '(A12, I8, A18)')   '(SE24.16e3, ', 1, '(", ",SE24.16e3))'
  write(format_V2, '(A12, I8, A18)')   '(SE24.16e3, ', 2*Nm_o+1, '(", ",SE24.16e3))'
  write(format_V3, '(A12, I8, A18)')   '(SE24.16e3, ', 2, '(", ",SE24.16e3))'
  write(format_V4, '(A12, I8, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'


write(*,*) shape(Ny1)
!  write(list_file, '(A)') 'f_end.dat'             !f
!  open(unit=709,file=list_file)
!  do i4 = 10:10
!  A_excit = (1+dble(i4)*0.1)*1d5*(Eg/3d0/hbar)/A_freq_para

!-------------------eigen vectors --------------------
  A = -1d0*TRANSPOSE(coul_mat(1, :, :))
  do i2 = 1, Ny
    A(i2, i2) = A(i2, i2) + y(i2)*y(i2)
  end do 
  LWORK = -1
  CALL ZGEEV(YES, YES, Ny, A, LDA, W, VL, LDVL,&
             VR, LDVR, WORK, LWORK, RWORK, INFO )
  LWORK = min( LWMAX, int( WORK( 1 ) ) )
  CALL ZGEEV(YES, YES, Ny, A, LDA, W, VL, LDVL,&
             VR, LDVR, WORK, LWORK, RWORK, INFO )
write(*,*) 'LWORK = ', LWORK
  A1 = -1d0*TRANSPOSE(coul_mat(2, :, :))
  LWORK = -1
  do i2 = 1, Ny
    A1(i2, i2) = A1(i2, i2) + y(i2)*y(i2)
  end do 
  CALL ZGEEV(YES, YES, Ny, A1, LDA, W1, VL1, LDVL,&
             VR1, LDVR, WORK, LWORK, RWORK, INFO )
  LWORK = min( LWMAX, int( WORK( 1 ) ) )
  CALL ZGEEV(YES, YES, Ny, A1, LDA, W1, VL1, LDVL,&
             VR1, LDVR, WORK, LWORK, RWORK, INFO )

  A2 = -1d0*TRANSPOSE(coul_mat(3, :, :))
  LWORK = -1
  do i2 = 1, Ny
    A2(i2, i2) = A2(i2, i2) + y(i2)*y(i2)
  end do 
  CALL ZGEEV(YES, YES, Ny, A2, LDA, W2, VL2, LDVL,&
             VR2, LDVR, WORK, LWORK, RWORK, INFO )
  LWORK = min( LWMAX, int( WORK( 1 ) ) )
  CALL ZGEEV(YES, YES, Ny, A2, LDA, W2, VL2, LDVL,&
             VR2, LDVR, WORK, LWORK, RWORK, INFO )

  A3 = -1d0*TRANSPOSE(coul_mat(4, :, :))
  LWORK = -1
  do i2 = 1, Ny
    A3(i2, i2) = A3(i2, i2) + y(i2)*y(i2)
  end do 
  CALL ZGEEV(YES, YES, Ny, A3, LDA, W3, VL3, LDVL,&
             VR3, LDVR, WORK, LWORK, RWORK, INFO )
  LWORK = min( LWMAX, int( WORK( 1 ) ) )
  CALL ZGEEV(YES, YES, Ny, A3, LDA, W3, VL3, LDVL,&
             VR3, LDVR, WORK, LWORK, RWORK, INFO )

!---------------end eigen vectors -------------------- 

!---------------sorting-------------------------------
  Etemp = real(W)
  Emax = maxval(Etemp)
  do Ndo = 1, Ny
    i2 = minloc(Etemp, 1)
    Etemp1(Ndo) = minval(Etemp, 1)
    VR_temp(:,Ndo) = VR(:,i2)
    VL_temp(:,Ndo) = VL(:,i2)
    Etemp(i2) = Emax+1d0
  end do
  W = Etemp1
  VR = VR_temp
  VL = VL_temp

  Etemp = real(W1)
  Emax = maxval(Etemp)
  do Ndo = 1, Ny
    i2 = minloc(Etemp, 1)
    Etemp1(Ndo) = minval(Etemp, 1)
    VR_temp(:,Ndo) = VR1(:,i2)
    VL_temp(:,Ndo) = VL1(:,i2)
    Etemp(i2) = Emax+1d0
  end do
  W1 = Etemp1
  VR1 = VR_temp
  VL1 = VL_temp

  Etemp = real(W2)
  Emax = maxval(Etemp)
  do Ndo = 1, Ny
    i2 = minloc(Etemp, 1)
    Etemp1(Ndo) = minval(Etemp, 1)
    VR_temp(:,Ndo) = VR2(:,i2)
    VL_temp(:,Ndo) = VL2(:,i2)
    Etemp(i2) = Emax+1d0
  end do
  W2 = Etemp1
  VR2 = VR_temp
  VL2 = VL_temp

  Etemp = real(W3)
  Emax = maxval(Etemp)
  do Ndo = 1, Ny
    i2 = minloc(Etemp, 1)
    Etemp1(Ndo) = minval(Etemp, 1)
    VR_temp(:,Ndo) = VR3(:,i2)
    VL_temp(:,Ndo) = VL3(:,i2)
    Etemp(i2) = Emax+1d0
  end do
  W3 = Etemp1
  VR3 = VR_temp
  VL3 = VL_temp

!-----------end sorting-------------------------------
!---------------Normalization-------------------------
  do Ndo = 1, Ny
    VL(:, Ndo) = VL(:, Ndo)/sqrt(real(sum(dy1*y*VL(:, Ndo)*conjg(VL(:, Ndo)))))
    VL1(:, Ndo) = VL1(:, Ndo)/sqrt(real(sum(dy1*y*VL1(:, Ndo)*conjg(VL1(:, Ndo)))))
    VL2(:, Ndo) = VL2(:, Ndo)/sqrt(real(sum(dy1*y*VL2(:, Ndo)*conjg(VL2(:, Ndo)))))
  end do
  do Ndo = 1, Ny
    VR(:, Ndo) = VR(:, Ndo)/sqrt(real(sum(dy1*y*VR(:, Ndo)*conjg(VR(:, Ndo)))))
    VR1(:, Ndo) = VR1(:, Ndo)/sqrt(real(sum(dy1*y*VR1(:, Ndo)*conjg(VR1(:, Ndo)))))
    VR2(:, Ndo) = VR2(:, Ndo)/sqrt(real(sum(dy1*y*VR2(:, Ndo)*conjg(VR2(:, Ndo)))))
  end do
  write(list_file, '(A)') 'oscillator_strength.dat'       !eigen values
  open(unit=700,file=list_file)
    write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
    do Ndo = 1, Ny
      osc_str(Ndo) = sum(y*dipole*VR(:, Ndo)*dy1/2d0/pi)
    end do
    osc_str = osc_str*conjg(osc_str)
    write(700, format_V) abs(osc_str)
  close(700)
!  VL1 = VR1
!  VL =VR
!  VL2 =VR2




  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  write(list_file, '(A)') 'p_r.dat'             !p
  open(unit=701,file=list_file)
  write(list_file2, '(A)') 'Et.dat'               !E(t)
  open(unit=702,file=list_file2)
  write(list_file, '(A)') 'pt_r.dat'           !p(t)
  open(unit=703,file=list_file)
  write(list_file, '(A)') 'pt_i.dat'           !p(t)
  open(unit=809,file=list_file)
  write(list_file, '(A)') 'time.dat'           !p(t)
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
  Ndo = 0.0d0
  do Ndo_in = 1, 2*Nm_o+1
    ft(Ndo_in) = sum(dy1*y*f(Ndo_in, :))   !calculate macroscopic density
    pt(Ndo_in) = sum(dy1*y*p(Ndo_in, :))/(2.0d0*pi)/a_B/a_B   !calculate macroscopic polarization
  end do
  J_THZ_t = sum(dy1*y*y*(f(Nm_o, :)+ f(Nm_o+2, :)))/(4.0d0*pi)   !calculate macroscopic density

  write(703, format_V2) real(pt)
  write(809, format_V2) aimag(pt)
  write(704, format_V2) abs(ft)
  write(711, format_V3) (Ndo*dt)+tstart_A    
  write(702, format_V1) real(Etime(dble(Ndo)))*dipole(1)
  write(705, format_V3) real(Atime(dble(Ndo))), aimag(Atime(dble(Ndo)))
  write(706, format_V3) real(J_THZ_t), aimag(J_THZ_t)

  do Ndo = 1, Nt
    
    !calculate p(0)

    call RK(E_freq, p_freq, p1_freq, f, p, decay, dble(Ndo), A_freq, J_THZ_freq, &              !input
            E_freq, p_freq, p1_freq, f, p, decay, A_freq, J_THZ_freq)                           !output
    do Ndo_in = 1, 2*Nm_o+1
!      ft(Ndo_in) = sum(dy1*y*f(Ndo_in, :))/(2.0d0*pi)   !calculate macroscopic density
      ft(Ndo_in) = sum(dy1*y*f(Ndo_in, :))*2.0d0   !calculate macroscopic density
      pt(Ndo_in) = sum(dy1*y*p(Ndo_in, :))/(2.0d0*pi)/a_B/a_B   !calculate macroscopic polarization
    end do
    J_THZ_t = sum(dy1*y*y*(f(Nm_o, :)+ f(Nm_o+2, :)))/(4.0d0*pi)   !calculate macroscopic density
    if(Ndo == i2*FLOOR(dble(Nt)/dble(Nt_RWA))) then
      write(703, format_V2) real(pt)
      write(809, format_V2) aimag(pt)
      write(704, format_V2) abs(ft)
      write(711, format_V3) (Ndo*dt)+tstart_A    

     
      write(702, format_V1) real(Etime(dble(Ndo)))*dipole(1)
      write(705, format_V3) real(Atime(dble(Ndo))), aimag(Atime(dble(Ndo)))
      write(706, format_V3) real(J_THZ_t), aimag(J_THZ_t)
!      p_proj1 = matmul(p(Nm_o+2, :)*y, conjg(VL1))*dy1/2d0/pi
!      p_proj = matmul(p(Nm_o+1, :)*y, conjg(VL))*dy1/2d0/pi
!      p_proj2 = matmul(p(Nm_o-1, :)*y, conjg(VL2))*dy1/2d0/pi
!      write(800, format_V4) real(p_proj)
!      write(800, format_V4) aimag(p_proj)
!      write(801, format_V4) real(p_proj1)
!      write(801, format_V4) aimag(p_proj1)
!      write(802, format_V4) real(p_proj2)
!      write(802, format_V4) aimag(p_proj2)
      i2 = i2+1
    end if
!    if(Ndo == i3*FLOOR(dble(Nt)/dble(Nt_RWA))) then
!      do i1 = 1, 2*Nm_o+1
!        write(701, format_V) real(p(i1, :))
!        write(710, format_V) aimag(p(i1, :))
!        write(707, format_V) real(f(i1, :))
!        write(708, format_V) aimag(f(i1, :))
!        write(807, format_V) real(decay(i1, :))
!        write(808, format_V) aimag(decay(i1, :))
!      end do
!        write(707, format_V) real(f(Nm_o+1, :))
!        write(708, format_V) aimag(f(Nm_o+1, :))

!      i3 = i3+1
!    end if
  end do

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
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Nm_o+1
    do i2 = 1, Ny
      write(700, format_V) real(coul_mat(i1, i2, :))
    end do
  end do
  close(700)
  write(list_file, '(A)') 'p_freq.dat'           !p_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) p_freq(i1)
  end do
  close(700)

  write(list_file, '(A)') 'p1_freq.dat'           !p_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) p1_freq(i1)
  end do
  close(700)
    
  write(list_file, '(A)') 'E_freq.dat'           !E_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) E_freq(i1)
  end do
  close(700)

  write(list_file, '(A)') 'A_freq.dat'           !E_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) A_freq(i1)
  end do
  close(700)

  write(list_file, '(A)') 'J_freq.dat'           !E_freq
  open(unit=700,file=list_file)
  write(format_V, '(A12, I6, A18)')   '(SE24.16e3, ', N_freq, '(", ",SE24.16e3))'
  do i1 = 1, N_freq
    write(700, format_V) J_THZ_freq(i1)
  end do
  close(700)

  write(list_file, '(A)') 'y.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) y
  close(700)

  write(list_file, '(A)') 'dy1.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) dy1
  close(700)

  write(list_file, '(A)') 'dipole.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) dipole
  close(700)


  write(list_file, '(A)') 'para_mat.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', 1, '(", ",SE24.16e3))'  
    write(700, format_V) t_end		!1
    write(700, format_V) dt*1d6   		!2
    write(700, *) Nt		!3
    write(700, format_V) ymax		!4
    write(700, *) Ny		!5
    write(700, format_V) A_excit*A_freq_para*1d-8		!6
    write(700, *) gamma		!7
    write(700, *) 2d0*pi/A_freq_para*1000d0		!8
    write(700, *) Nt_RWA		!9
    write(700, *) tstart_A		!10
    write(700, *) sigmat_A		!11
    write(700, *) Nm_o		!12
    write(700, *) N_freq		!13
    write(700, *) delta_1s		!14
    write(700, *) decay_m		!15
    write(700, *) gamma		!16
    write(700, *) E_B		!17
    write(700, *) a_B		!18
    write(700, *) Eg		!19
    write(700, *) N_fine		!20
    write(700, *) Nphi		!21
  close(700)
  write(list_file, '(A)') 'eigenvalues_s.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) real(W)
  close(700)
  write(list_file, '(A)') 'eigenvalues_p.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) real(W1)
  close(700)
  write(list_file, '(A)') 'eigenvalues_d.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) real(W2)
  close(700)
  write(list_file, '(A)') 'eigenvalues_f.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) real(W3)
  close(700)

  write(list_file, '(A)') 'test.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) test
  close(700)
  write(list_file, '(A)') 'eigenvec_1s.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) real(VR(:,2))
  close(700)
  write(list_file, '(A)') 'eigenvec_2p.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
      write(700, format_V) real(VR1(:,1))
  close(700)
  write(list_file, '(A)') 'eigenvec_s.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Ny
      write(700, format_V) real(VR(:,i1))
  end do
  close(700)
  write(list_file, '(A)') 'eigenvec_p.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Ny
      write(700, format_V) real(VR1(:,i1))
  end do
  close(700)
  write(list_file, '(A)') 'eigenvec_d.dat'       !eigen values
  open(unit=700,file=list_file)
  write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
  do i1 = 1, Ny
      write(700, format_V) real(VR2(:,i1))
  end do
  close(700)
!----------------end output------------------
  call output_params
end program main
