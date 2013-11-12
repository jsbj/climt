<<<<<<< HEAD
subroutine driver(jm, im, srf_rad_flx, srf_lat_flx, srf_sen_flx, Hslab, rowl, Cl, dt, &
=======
subroutine driver(jm, im, srf_rad_flx, srf_lat_flx, srf_sen_flx, qflx, Hslab, rowl, Cl, dt, &
>>>>>>> 1055422fefa06f7092b71f9ab4c7a66a6ba25840
                             Tsinc, Tsdot )

implicit none

! Input
integer, intent(in)                   :: jm,im
<<<<<<< HEAD
real(8), intent(in), dimension(jm,im) :: srf_rad_flx, srf_lat_flx, srf_sen_flx
=======
real(8), intent(in), dimension(jm,im) :: srf_rad_flx, srf_lat_flx, srf_sen_flx, qflx
>>>>>>> 1055422fefa06f7092b71f9ab4c7a66a6ba25840
real(8), intent(in)                   :: Hslab, rowl, Cl, dt
!f2py intent(in,hide)  jm,im

! Output
real(8), intent(out), dimension(jm,im) :: Tsinc, Tsdot

! Local
real(8) :: invCslab

invCslab = 1. / (Hslab*rowl*Cl)
<<<<<<< HEAD
Tsdot = (srf_rad_flx + srf_lat_flx + srf_sen_flx) * invCslab
=======
Tsdot = (srf_rad_flx + srf_lat_flx + srf_sen_flx + qflx) * invCslab
>>>>>>> 1055422fefa06f7092b71f9ab4c7a66a6ba25840
Tsinc = 2.*dt*Tsdot

end subroutine driver
