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
    real*16                                 ::a(Nphi), b(Nphi), c(Nphi), w11(Nphi)   
!a:integer series   b&w11: Gauss-Chebyshev grid
    integer                                          ::Ndo, Ndo_in, Ndo_m, m = 0   
!Ndo & Ndo_in:integer number for do loop, m:order of circular harmonics 
    real*16                                 ::fine_grid(N_fine)
!fine_grid: finer grid for singularity removal
    real*16                                 ::b_mid1, b_mid2(Nphi), b_mid3(Nphi)  
!help terms to help calculating Gauss-Chebyshev grid
    do Ndo = 1, Nphi
      a(Ndo) = real(Ndo,16)             
      !a:integer series from 1 to Nphi
    end do
    do Ndo = 1, N_fine
      fine_grid(Ndo) = (real(Ndo,16)-0.50q0)*dy_fine  
      !fine_grid: finer grid for singularity removal
    end do
    b_mid1 = ( real(Nphi,16) + 1.0q0 ) / (2.0q0 * pi )
    b_mid2 = sin ( (2.0q0 * pi * a ) / ( real(Nphi,16)+ 1.0q0 ))
    b_mid3 =  ( a - b_mid1 * b_mid2)
    b = pi / ( real(Nphi,16) + 1.0q0 ) * b_mid3
    w11 = pi / ( real(Nphi,16) + 1.0q0 ) * ( 1.0q0 - cos( 2.0q0 * pi * a / ( real(Nphi,16) + 1.0q0))) 
    !b&w11: Gauss-Chebyshev grid
    do Ndo_m = 1, Nm_o+1
      do Ndo = 1, Ny
        y_fine = y(Ndo) - dy/2.0q0 + fine_grid   
        !calculate the finer grid for diagonal terms
!------------calculate the Coulomb matrix without removal of singularity-------------
        do Ndo_in = 1, Nphi
          coul_mat(Ndo_m, :, Ndo)=coul_mat(Ndo_m, :, Ndo)+(1.0q0/(sqrt(y(Ndo)*y(Ndo)+y*y -2.0q0*y*y(Ndo)*cos(b(Ndo_in)))))*2.0q0*y*w11(Ndo_in)*cos(real(Ndo_m-1,16)*b(Ndo_in))
        end do
!-----------end calculation----------------------------------------------------------
        coul_mat(Ndo_m, Ndo, Ndo) = 0.0q0     
       !eliminate the diagonal terms and calculate them as followed
!-----------calculate the diagonal terms with removal of singularity---------------------
        do Ndo_in = 1, N_fine    
          coul_mat(Ndo_m, Ndo,Ndo)=coul_mat(Ndo_m, Ndo, Ndo) + sum( 1.0q0/ sqrt( y(Ndo) *y(Ndo) + y_fine(Ndo_in) *y_fine(Ndo_in) - 2.0q0 * y_fine(Ndo_in) * y(Ndo) * cos( b ) ) / real(N_fine,16) * 2.0q0 *y_fine(Ndo_in) * w11* cos( real(Ndo_m-1,16)*b ))
        end do  
!-----------end calculation-----------------------------------------------------------
      end do
    end do
    coul_mat = coul_mat/pi*dy   
    !multiply with constants that missed in the do loop
    coul_mat = 0q0
  end subroutine coul_matrix

end module coul_mat_module
