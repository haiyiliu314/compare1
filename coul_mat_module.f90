!08/21/2017 Revision
!Add m as the order of circular harmonics
module coul_mat_module
  use constants
  use N_of_grids
  implicit none
  contains
  subroutine coul_matrix
  !07/03/2017 creation
  !calculate Coulomb matrix(non-symmetric, i.e. already contains the integration term, i.e. 2pi*dy*y*V_{y-y'})
    implicit none
    integer, parameter                               ::N_read = 101
    double precision                                 ::a(Nphi), b(Nphi), c(Nphi), w11(Nphi)   
!a:integer series   b&w11: Gauss-Chebyshev grid
    integer                                          ::Ndo, Ndo_in, Ndo_in_in, Ndo_m, m = 0   
!Ndo & Ndo_in:integer number for do loop, m:order of circular harmonics 
    double precision                                 ::fine_grid(N_fine)
!fine_grid: finer grid for singularity removal
    double precision                                 ::b_mid1, b_mid2(Nphi), b_mid3(Nphi)  
!help terms to help calculating Gauss-Chebyshev grid
    double precision                                 ::y_para(N_read), ghh_para(N_read), gch_para(N_read), gcc_para(N_read), gcc_conf(Ny), gch_conf(Ny), ghh_conf(Ny), y_leng(Ny), y_leng1, gcc_conf1, gch_conf1, ghh_conf1
    integer                                          ::round_y(Ny), round_y1
    do Ndo = 1, Nphi
      a(Ndo) = dble(Ndo)             
      !a:integer series from 1 to Nphi
    end do
    do Ndo = 1, N_fine
      fine_grid(Ndo) = (dble(Ndo)-0.50d0)*dy_fine  
      !fine_grid: finer grid for singularity removal
    end do
  write(list_file, '(A)') 'pk2.dat'           !pk
  open(unit=700,file=list_file)
    open(unit = 100, file = 'y_para.txt', status = 'old', action = 'read')
    open(unit = 101, file = 'ghh_para.txt', status = 'old', action = 'read')
    open(unit = 102, file = 'gch_para.txt', status = 'old', action = 'read')
    open(unit = 103, file = 'gcc_para.txt', status = 'old', action = 'read')
    do Ndo = 1,N_read
      read(100,*) y_para(Ndo)
      read(101,*) ghh_para(Ndo)
      read(102,*) gch_para(Ndo)
      read(103,*) gcc_para(Ndo)
    end do


    b_mid1 = ( dble(Nphi) + 1.0d0 ) / (2.0d0 * pi )
    b_mid2 = dsin ( (2.0d0 * pi * a ) / ( dble(Nphi)+ 1.0d0 ))
    b_mid3 =  ( a - b_mid1 * b_mid2)
    b = pi / ( dble(Nphi) + 1.0d0 ) * b_mid3
    w11 = pi / ( dble(Nphi) + 1.0d0 ) * ( 1.0d0 - dcos( 2.0d0 * pi * a / ( dble(Nphi) + 1.0d0))) 
    !b&w11: Gauss-Chebyshev grid
    do Ndo_m = 1, Nm_o+1
      do Ndo = 1, Ny
        y_fine = y(Ndo) - dy/2.0d0 + fine_grid   

        !calculate the finer grid for diagonal terms
!------------calculate the Coulomb matrix without removal of singularity-------------
        do Ndo_in = 1, Nphi
          y_leng = dsqrt(y(Ndo)*y(Ndo)+y*y -2.0d0*y*y(Ndo)*dcos(b(Ndo_in)))
          round_y = int(y_leng*dble(N_read-1)/y_para(N_read)) !Nread: the number of elements in y_para.txt
          gcc_conf = 0.0d0
          gch_conf = 0.0d0
          ghh_conf = 0.0d0
          do Ndo_in_in = 1,Ny
            gcc_conf(Ndo_in_in) = gcc_para(round_y(Ndo_in_in)+1) + (y_leng(Ndo_in_in) - y_para(round_y(Ndo_in_in)+1))/(y_para(round_y(Ndo_in_in)+2)-y_para(round_y(Ndo_in_in)+1))*(gcc_para(round_y(Ndo_in_in)+2)-gcc_para(round_y(Ndo_in_in)+1))
            gch_conf(Ndo_in_in) = gch_para(round_y(Ndo_in_in)+1) + (y_leng(Ndo_in_in) - y_para(round_y(Ndo_in_in)+1))/(y_para(round_y(Ndo_in_in)+2)-y_para(round_y(Ndo_in_in)+1))*(gch_para(round_y(Ndo_in_in)+2)-gch_para(round_y(Ndo_in_in)+1))
            ghh_conf(Ndo_in_in) = ghh_para(round_y(Ndo_in_in)+1) + (y_leng(Ndo_in_in) - y_para(round_y(Ndo_in_in)+1))/(y_para(round_y(Ndo_in_in)+2)-y_para(round_y(Ndo_in_in)+1))*(ghh_para(round_y(Ndo_in_in)+2)-ghh_para(round_y(Ndo_in_in)+1))
          end do
          coul_mat(Ndo_m, :, Ndo)=coul_mat(Ndo_m, :, Ndo)+(1.0d0/y_leng)*2.0d0*y*w11(Ndo_in)*dcos(dble(Ndo_m-1)*b(Ndo_in)) /pi*dy   
          coul_mat_cc(Ndo_m, :, Ndo)=coul_mat_cc(Ndo_m, :, Ndo)+(1.0d0/y_leng)*2.0d0*y*w11(Ndo_in)*dcos(dble(Ndo_m-1)*b(Ndo_in))*gcc_conf/pi*dy   
          coul_mat_ch(Ndo_m, :, Ndo)=coul_mat_ch(Ndo_m, :, Ndo)+(1.0d0/y_leng)*2.0d0*y*w11(Ndo_in)*dcos(dble(Ndo_m-1)*b(Ndo_in))*gch_conf/pi*dy   
          coul_mat_hh(Ndo_m, :, Ndo)=coul_mat_hh(Ndo_m, :, Ndo)+(1.0d0/y_leng)*2.0d0*y*w11(Ndo_in)*dcos(dble(Ndo_m-1)*b(Ndo_in))*ghh_conf/pi*dy   
        end do

!-----------end calculation----------------------------------------------------------
        coul_mat(Ndo_m, Ndo, Ndo) = 0.0d0     
        coul_mat_cc(Ndo_m, Ndo, Ndo) = 0.0d0     
        coul_mat_ch(Ndo_m, Ndo, Ndo) = 0.0d0     
        coul_mat_hh(Ndo_m, Ndo, Ndo) = 0.0d0     
       !eliminate the diagonal terms and calculate them as followed
!-----------calculate the diagonal terms with removal of singularity---------------------
        do Ndo_in = 1, N_fine 
          gcc_conf1 = 0.0d0
          gch_conf1 = 0.0d0
          ghh_conf1 = 0.0d0
          do Ndo_in_in = 1,Nphi
            y_leng1 =  dsqrt( y(Ndo) *y(Ndo) + y_fine(Ndo_in) *y_fine(Ndo_in) - 2.0d0 * y_fine(Ndo_in) * y(Ndo) * dcos( b(Ndo_in_in) ) )
            round_y1 = int(y_leng1*N_read/y_para(N_read))
            gcc_conf1 = gcc_para(round_y1+1) + (y_leng1 - y_para(round_y1+1))/(y_para(round_y1+2)-y_para(round_y1+1))*(gcc_para(round_y1+2)-gcc_para(round_y1+1))
            gch_conf1 = gch_para(round_y1+1) + (y_leng1 - y_para(round_y1+1))/(y_para(round_y1+2)-y_para(round_y1+1))*(gch_para(round_y1+2)-gch_para(round_y1+1))
            ghh_conf1 = ghh_para(round_y1+1) + (y_leng1 - y_para(round_y1+1))/(y_para(round_y1+2)-y_para(round_y1+1))*(ghh_para(round_y1+2)-ghh_para(round_y1+1))
            coul_mat(Ndo_m, Ndo, Ndo)=coul_mat(Ndo_m, Ndo, Ndo) + (1.0d0/y_leng1 * 2.0d0 *y_fine(Ndo_in) * w11(Ndo_in_in)* dcos( dble(Ndo_m-1)*b(Ndo_in_in) ))/pi*dy/ dble(N_fine)  
            coul_mat_cc(Ndo_m, Ndo, Ndo)=coul_mat_cc(Ndo_m, Ndo, Ndo) + (1.0d0/y_leng1 * 2.0d0 *y_fine(Ndo_in) * w11(Ndo_in_in)* dcos( dble(Ndo_m-1)*b(Ndo_in_in) ))*gcc_conf1/pi*dy/ dble(N_fine)   
            coul_mat_ch(Ndo_m, Ndo, Ndo)=coul_mat_ch(Ndo_m, Ndo, Ndo) + (1.0d0/y_leng1 * 2.0d0 *y_fine(Ndo_in) * w11(Ndo_in_in)* dcos( dble(Ndo_m-1)*b(Ndo_in_in) ))*gch_conf1/pi*dy/ dble(N_fine)   
            coul_mat_hh(Ndo_m, Ndo, Ndo)=coul_mat_hh(Ndo_m, Ndo, Ndo) + (1.0d0/y_leng1 * 2.0d0 *y_fine(Ndo_in) * w11(Ndo_in_in)* dcos( dble(Ndo_m-1)*b(Ndo_in_in) ))*ghh_conf1/pi*dy/ dble(N_fine)   
!            write(format_V, '(A12, I4, A18)')   '(SE24.16e3, ', Ny, '(", ",SE24.16e3))'
!            write(700, format_V) y_leng1
          end do
        end do  
!-----------end calculation-----------------------------------------------------------
      end do
    end do
  close(700)
!    coul_mat = coul_mat/pi*dy   
!    coul_mat_cc = coul_mat_cc/pi*dy   
!    coul_mat_ch = coul_mat_ch/pi*dy   
!    coul_mat_hh = coul_mat_hh/pi*dy   
    !multiply with constants that missed in the do loop
!    coul_mat = 0d0
!    coul_mat = coul_mat*10d0

  end subroutine coul_matrix

end module coul_mat_module
