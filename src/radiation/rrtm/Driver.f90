! see _rrtm_radiation for the python that prepares these arguments...
subroutine driver
    (ncol, nlay, icld, idrv, play, plev, tlay, tlev, tsfc, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, &
    o2vmr, cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, aldif, aldir, asdif, asdir, emis, adjes, coszen, scon, &
    inflgsw, iceflgsw, liqflgsw, cldfmcl_sw, taucmcl_sw, ssacmcl_sw, asmcmcl_sw, fsfcmcl_sw, ciwpmcl_sw, clwpmcl_sw, &
    reicmcl_sw, relqmcl_sw, inflglw, iceflgslw, liqflglw, cldfmcl_lw, ciwpmcl_lw, clwpmcl_lw, reicmcl_lw, &
    relqmcl_lw, taucmcl_lw, tauaer_sw, ssaaer_sw, asmaer_sw, ecaer_sw, tauaer_lw, swuflx, swdflx, swhr, swuflxc, &
    swdflxc, swhrc, uflx, dflx, hr, uflxc, dflxc, hrc)
! Modules
    use rrtmg_lw_rad, only: rrtmg_lw
    use rrtmg_sw_rad, only: rrtmg_sw
    use parkind, only : im => kind_im, rb => kind_rb
    use parrrtm, only : nbndlw, ngptlw, maxxsec, mxmol
    use parrrsw, only : nbndsw, ngptsw, naerec, nstr, nmol, mxmol, jpband, jpb1, jpb2
    
! Input
    integer(kind=im), intent(in) :: ncol
    integer(kind=im), intent(in) :: nlay
    integer(kind=im), intent(in) :: icld
    integer(kind=im), intent(in) :: idrv
    real(kind=rb), intent(in) :: play(ncol,nlay)
    real(kind=rb), intent(in) :: plev(ncol,nlay+1)
    real(kind=rb), intent(in) :: tlay(ncol,nlay)
    real(kind=rb), intent(in) :: tlev(ncol,nlay+1)
    real(kind=rb), intent(in) :: tsfc(ncol)
    real(kind=rb), intent(in) :: h2ovmr(ncol,nlay)
    real(kind=rb), intent(in) :: o3vmr(ncol,nlay)
    real(kind=rb), intent(in) :: co2vmr(ncol,nlay)
    real(kind=rb), intent(in) :: ch4vmr(ncol,nlay)
    real(kind=rb), intent(in) :: n2ovmr(ncol,nlay)
    real(kind=rb), intent(in) :: o2vmr(ncol,nlay)
    real(kind=rb), intent(in) :: cfc11vmr(ncol,nlay)
    real(kind=rb), intent(in) :: cfc12vmr(ncol,nlay)
    real(kind=rb), intent(in) :: cfc22vmr(ncol,nlay)
    real(kind=rb), intent(in) :: ccl4vmr(ncol,nlay)
    real(kind=rb), intent(in) :: aldif(ncol)
    real(kind=rb), intent(in) :: aldir(ncol)
    real(kind=rb), intent(in) :: asdif(ncol)
    real(kind=rb), intent(in) :: asdir(ncol)
    real(kind=rb), intent(in) :: emis(ncol,nbndlw)
    real(kind=rb), intent(in) :: adjes
    real(kind=rb), intent(in) :: coszen(ncol)
    real(kind=rb), intent(in) :: scon
    integer(kind=im), intent(in) :: inflgsw
    integer(kind=im), intent(in) :: iceflgsw
    integer(kind=im), intent(in) :: liqflgsw
    real(kind=rb), intent(in) :: cldfmcl_sw(ngptsw,ncol,nlay)
    real(kind=rb), intent(in) :: taucmcl_sw(ngptsw,ncol,nlay)
    real(kind=rb), intent(in) :: ssacmcl_sw(ngptsw,ncol,nlay)
    real(kind=rb), intent(in) :: asmcmcl_sw(ngptsw,ncol,nlay)
    real(kind=rb), intent(in) :: fsfcmcl_sw(ngptsw,ncol,nlay)
    real(kind=rb), intent(in) :: ciwpmcl_sw(ngptsw,ncol,nlay)
    real(kind=rb), intent(in) :: clwpmcl_sw(ngptsw,ncol,nlay)
    real(kind=rb), intent(in) :: reicmcl_sw(ncol,nlay)
    real(kind=rb), intent(in) :: relqmcl_sw(ncol,nlay)
    integer(kind=im), intent(in) :: inflglw
    integer(kind=im), intent(in) :: iceflgslw
    integer(kind=im), intent(in) :: liqflglw
    real(kind=rb), intent(in) :: cldfmcl_lw(ngptlw,ncol,nlay)
    real(kind=rb), intent(in) :: ciwpmcl_lw(ngptlw,ncol,nlay)
    real(kind=rb), intent(in) :: clwpmcl_lw(ngptlw,ncol,nlay)
    real(kind=rb), intent(in) :: reicmcl_lw(ncol,nlay)
    real(kind=rb), intent(in) :: relqmcl_lw(ncol,nlay)
    real(kind=rb), intent(in) :: taucmcl_lw(ngptlw,ncol,nlay)
    real(kind=rb), intent(in) :: tauaer_sw(ncol,nlay,nbndsw)
    real(kind=rb), intent(in) :: ssaaer_sw(ncol,nlay,nbndsw)
    real(kind=rb), intent(in) :: asmaer_sw(ncol,nlay,nbndsw)
    real(kind=rb), intent(in) :: ecaer_sw(ncol,nlay,nbndsw)    
    real(kind=rb), intent(in) :: tauaer_lw(ncol,nlay,nbndlw)
    ! Output

    ! SW
    real(kind=rb), intent(out) :: swuflx(ncol,nlay+1)       ! Total sky shortwave upward flux (W/m2)
    real(kind=rb), intent(out) :: swdflx(ncol,nlay+1)       ! Total sky shortwave downward flux (W/m2)
    real(kind=rb), intent(out) :: swhr(ncol,nlay)         ! Total sky shortwave radiative heating rate (K/d)
    real(kind=rb), intent(out) :: swuflxc(ncol,nlay+1)      ! Clear sky shortwave upward flux (W/m2)
    real(kind=rb), intent(out) :: swdflxc(ncol,nlay+1)      ! Clear sky shortwave downward flux (W/m2)
    real(kind=rb), intent(out) :: swhrc(ncol,nlay)        ! Clear sky shortwave radiative heating rate (K/d)
                                                
    ! LW
    real(kind=rb), intent(out) :: uflx(ncol,nlay+1)         ! Total sky longwave upward flux (W/m2)
    real(kind=rb), intent(out) :: dflx(ncol,nlay+1)         ! Total sky longwave downward flux (W/m2)
    real(kind=rb), intent(out) :: hr(ncol,nlay)           ! Total sky longwave radiative heating rate (K/d)
    real(kind=rb), intent(out) :: uflxc(ncol,nlay+1)        ! Clear sky longwave upward flux (W/m2)
    real(kind=rb), intent(out) :: dflxc(ncol,nlay+1)        ! Clear sky longwave downward flux (W/m2)
    real(kind=rb), intent(out) :: hrc(ncol,nlay)          ! Clear sky longwave radiative heating rate (K/d)
                                                
end subroutine driver
