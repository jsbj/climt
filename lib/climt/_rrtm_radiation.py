from scipy.interpolate import interp1d
from math import cos
import _rrtm_radiation_fortran

# wavenumber bands used by RRTM:
SW_BANDS = range(14)
LW_BANDS = range(16)


INPUTS = [
            # 'do_sw', #     0  Shortwave switch     (integer)       1   1 / 0 => do / do not compute SW        
            # 'do_lw', #     0  Longwave switch      (integer)       1   1 / 0 => do / do not compute LW       
            'p', #       1-3  Atmospheric pressure     mb               Default is equispaced 0-ps. p[0] is top level
            'T', #       1-3  Temperature              K        283.15  Isothermal
            'ps', #      0-2  Surface pressure         mb       1000.
            'Ts', #      0-2  Surface temperature      K        283.15
            'q', #       1-3  Specific humidity        g/kg     1.e-5
            'o3', #      1-3  Ozone mass mix. rat.     kg/kg            Default obtained by interpolating a tropical data profile
            'co2', #       0  CO2                  ppmv          330.
            'ch4', #       0  CH4                  ppmv            0.
            'n2o', #       0  N2O                  ppmv            0.
            'o2', #       0  O2                    volume mixing ratio
            'cfc11', #     0  CFC11                ppmv            0.
            'cfc12', #     0  CFC12                ppmv            0.
            'cfc22', #     0  CFC22                ppmv            0.
            'ccl4',  #        CCl4                volume mixing ratio
            'aldif', #   0-2  Diffuse near-IR (SW) albedo        (frac)   0.07
            'aldir', #   0-2  Direct near-IR (SW) albedo         (frac)   0.07
            'asdif', #   0-2  Diffuse UV+vis alb       (frac)   0.07
            'asdir', #   0-2  Direct UV+vis alb       (frac)   0.07
            'lw_surface_emissivity', # should have len(LW_BANDS) members...see above
            'zen', #     0-2  Solar zenith angle       dgr      72.2    Daily-mean on equator at equinox
            # 'calday', #    0  Calendar day         (float)      80.5    Insolation computed at specified
            # 'orb_yr', #    0  Orbital year         (integer)    1995   Year used to compute orbital params
            # 'avg', #       0  Insolation average   (string)   'daily'  Choices are: 'inst', 'daily', 'annual'
            # 'lat', #     0-1  Latitude                 dgr      0.         day and lat/lon if  solin 
            # 'lon', #     0-1  Longitude                dgr      0.         and zen are NOT specified
            # 'solin', #   0-2  Insolation               W/m2     417.4   Daily-mean on equator at equinox
            'scon', #      0  Solar constant       W m-2        1367.
            
            # 'tauvis', #    0  Aerosol opt. depth   (float)         0.   CCM3 only
            # 'tau_inf', #   0  Total opt. depth        -            1.   Greygas scheme only
            # 'alpha_greygas', # 0  Tau shape parameter   -          1.   Greygas scheme only  
            'cldf', #    1-3  Cloud fraction           frac     0.:
            # 'in_cld', #    0  Cloud water path flag     -       0       0 / 1 => grid avg / in-cloud water paths (CAM3 only)
            'cloud_single_scattering_albedo',
            'cloud_asymmetry_parameter',
            'cloud_forward_scattering_fraction',
            'r_liq', #   1-3  Drop radius, liquid      micron   10.   
            'r_ice', #   1-3  Drop radius, ice         micron   30.   
            'clwp', #    1-3  Cloud liquid water path  g/m2     0.     
            'ciwp', #    1-3  Cloud ice water path     g/m2     -99.   If not passed explicitly, ice frac computed internally (CAM3 only)
            # 'flus' #    1-3  Upwelling LW at surface  W/m2     -99.   If not passed explicitly, computed from Ts using emiss=1 (CAM3 only)
            'tauaer_sw', # Aerosol optical depth (iaer=10 only), Dimensions: (ncol,nlay,nbndsw), (non-delta scaled)
            'ssaaer_sw', # Aerosol single scattering albedo (iaer=10 only), Dimensions: (ncol,nlay,nbndsw), (non-delta scaled)      
            'asmaer_sw', # Aerosol asymmetry parameter (iaer=10 only), Dimensions: (ncol,nlay,nbndsw), (non-delta scaled)
            'tauaer_lw' # Aerosol optical depth (iaer=10 only), Dimensions: (ncol,nlay,nbndlw), (non-delta scaled)
]

def driver(*args):
    # gotta translate between the APIs:
    climt_inputs = dict(zip(INPUTS, args))

    number_of_layers = len(climt_inputs['T'])
    
    
    climt_inputs['p'] = [p[0][0] for p in climt_inputs['p']]
    climt_inputs['T'] = [T[0][0] for T in climt_inputs['T']]
    
    if 'pbound' not in climt_inputs:
        interpolated_p = interp1d(range(number_of_layers), climt_inputs['p'])
    if 'Tbound' not in climt_inputs:
        interpolated_T = interp1d(range(number_of_layers), climt_inputs['T'])
    
    clouds = 1 if 'cldf' in climt_inputs else 0
    
    for key in ['co2', 'ch4', 'n2o', 'o2', 'cfc11', 'cfc12', 'cfc22', 'ccl4']:
        if not hasattr(climt_inputs[key], '__iter__'):
            climt_inputs[key] = [climt_inputs[key]] * number_of_layers
    for key in ['lw_surface_emissivity']:
        if not hasattr(climt_inputs[key], '__iter__'):
            climt_inputs[key] = [climt_inputs[key]] * len(LW_BANDS)

    rrtm_inputs = [
        # GENERAL, used in both SW and LW
        'ncol': 1, # number of columns
        'nlay': number_of_layers,
        'icld': clouds, # Cloud overlap method, 0: Clear only, 1: Random, 2: Maximum/random, 3: Maximum
        'idrv': 0, # whether to also calculate the derivative of flux with respect to surface temp
        'play': climt_inputs['p'], # pressure in each layer
        'plev': [climt_inputs['ps']] + [interpolated_p(i + .5) for i in range(number_of_layers - 1)] + [2 * climt_inputs['p'][-1] - interpolated_p(number_of_layers - 1.5)], # pressure at boundaries of each layer, with a linear extrapolation for the "top"
        'tlay': climt_inputs['T'], # temperature in each layer
        'tlev': [climt_inputs['Ts']] + [interpolated_p(i + .5) for i in range(number_of_layers - 1)] + [2 * climt_inputs['T'][-1] - interpolated_p(number_of_layers - 1.5)], # temperature at boundaries of each layer, with a linear extrapolation for the "top"
        'tsfc': climt_inputs['Ts'],
        # GASES, used in both SW and LW
        'h2ovmr': [((q/1000.)/(1. - (q/1000.)))*1.607793 for q in climt_inputs['q']], # convert from g/kg to volume mixing ration using molecular weight of dry air / water vapor
        'o3vmr': [o3 * 0.603428 for o3 in climt_inputs['o3']], # convert from kg/kg to volume mixing ratio using molecular weight of dry air / ozone
        'co2vmr': [co2 / 1.e6 for co2 in climt_inputs['co2']],
        'ch4vmr': [ch4 / 1.e6 for ch4 in climt_inputs['ch4']],
        'n2ovmr': [n2o / 1.e6 for n2o in climt_inputs['n2o']],
        'o2vmr': climt_inputs['o2'],
        'cfc11vmr': [cfc11 / 1.e6 for cfc11 in climt_inputs['cfc11']],
        'cfc12vmr': [cfc12 / 1.e6 for cfc12 in climt_inputs['cfc12']],
        'cfc22vmr': [cfc22 / 1.e6 for cfc22 in climt_inputs['cfc22']],
        'ccl4vmr': climt_inputs['ccl4'],
        # SURFACE OPTICAL PROPERTIES
        # SW
        'aldif': climt_inputs['aldif'],
        'aldir': climt_inputs['aldir'],
        'asdif': climt_inputs['asdif'],
        'asdir': climt_inputs['asdir'],
        # LW
        'emis': [1 - emis for emis in climt_inputs['lw_surface_emissivity']],
        # THE SUN - SW
        # 'dyofyr': # day of the year, used to get Earth/Sun distance (if not adjes)
        'adjes': 1., # flux adjustment for earth/sun distance (if not dyofyr)
        'coszen': cos(climt_inputs['zen']), # cosine of the solar zenith angle
        'scon': climt_inputs['scon'], # solar constant
        # CLOUDS, SW see http://www.arm.gov/publications/proceedings/conf16/extended_abs/iacono_mj.pdf
        'inflgsw': 2, # Flag for cloud optical properties
            # INFLAG = 0 direct specification of optical depths of clouds;
            #            cloud fraction and cloud optical depth (gray) are
            #            input for each cloudy layer
            #        = 1 calculation of combined ice and liquid cloud optical depths (gray) 
            #            as in CCM2; cloud fraction and cloud water path are input for
            #            each cloudy layer.  
            #        = 2 calculation of separate ice and liquid cloud optical depths, with
            #            parameterizations determined by values of ICEFLAG and LIQFLAG. 
            #            Cloud fraction, cloud water path, cloud ice fraction, and
            #            effective ice radius are input for each cloudy layer for all 
            #            parameterizations.  If LIQFLAG = 1, effective liquid droplet radius
            #            is also needed. 
        'iceflgsw': 0, # Flag for ice particle specification
            #             ICEFLAG = 0 the optical depths (gray) due to ice clouds are computed as in CCM3.
            #                     = 1 the optical depths (non-gray) due to ice clouds are computed as closely as
            #                         possible to the method in E.E. Ebert and J.A. Curry, JGR, 97, 3831-3836 (1992).
            #                     = 2 the optical depths (non-gray) due to ice clouds are computed by a method
            #                         based on the parameterization used in the radiative transfer model Streamer
            #                         (reference: J. Key, Streamer User's Guide, Technical Report 96-01, Boston
            #                         University, 85 pp. (1996)), which is closely related to the parameterization
            #                         of water clouds due to Hu and Stamnes (see below).
            #             = 3 the optical depths (non-gray) due to ice clouds are computed by a method
            # based on the parameterization given in Fu et al., J. Clim.,11,2223-2237 (1998).
        'liqflgsw': 1, # Flag for liquid droplet specification
            # LIQFLAG = 0 the optical depths (gray) due to water clouds are computed as in CCM3.
            #         = 1 the optical depths (non-gray) due to water clouds are computed by a method
            #             based on the parameterization of water clouds due to Y.X. Hu and K. Stamnes,
            #             J. Clim., 6, 728-742 (1993).
        'cldfmcl_sw': [[climt_inputs['cldf']]] * 112, # Cloud fraction; 112 is the "g-interval" for sw, ngptsw, set in parrrsw.f90
        'taucmcl_sw': [[None]] * 112, # In-cloud optical depth [IS THIS ONE NEEDED GIVEN THE OTHERS?]
        'ssacmcl_sw': [[climt_inputs['cloud_single_scattering_albedo']]] * 112, # In-cloud single scattering albedo
        'asmcmcl_sw': [[climt_inputs['cloud_asymmetry_parameter']]] * 112, # In-cloud asymmetry parameter
        'fsfcmcl_sw': [[climt_inputs['cloud_forward_scattering_fraction']]] * 112, # In-cloud forward scattering fraction (delta function pointing forward "forward peaked scattering")
        'ciwpmcl_sw': [[climt_inputs['ciwp']]] * 112, # in-cloud ice water path (g/m2)
        'clwpmcl_sw': [[climt_inputs['clwp']]] * 112, # in-cloud liquid water path (g/m2)
        'reicmcl_sw': [climt_inputs['r_ice']], # Cloud ice particle effective size (microns)
        'relqmcl_sw': [climt_inputs['r_liq']], # Cloud water drop effective radius (microns)
        'inflglw': 2,
        'iceflgslw': 0,
        'liqflglw': 1,
        'cldfmcl_lw': [[climt_inputs['cldf']]] * 140, # Cloud fraction; 140 is the "g-interval" for lw, ngptlw, set in parrrtm.f90
        'ciwpmcl_lw': [[climt_inputs['ciwp']]] * 140, # in-cloud ice water path (g/m2)
        'clwpmcl_lw': [[climt_inputs['clwp']]] * 140, # in-cloud liquid water path (g/m2)        
        'reicmcl_lw': [climt_inputs['r_ice']],    #  Cloud ice particle effective size (microns)
                                                              # specific definition of reicmcl depends on setting of iceflglw:
                                                              # iceflglw = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                              #               r_ec must be >= 10.0 microns
                                                              # iceflglw = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                              #               r_ec range is limited to 13.0 to 130.0 microns
                                                              # iceflglw = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                              #               r_k range is limited to 5.0 to 131.0 microns
                                                              # iceflglw = 3: generalized effective size, dge, (Fu, 1996),
                                                              #               dge range is limited to 5.0 to 140.0 microns
                                                              #               [dge = 1.0315 * r_ec]
        'relqmcl_lw': [climt_inputs['r_liq']], # Cloud water drop effective radius (microns)
        'taucmcl_lw': [[None]] * 140, # In-cloud optical depth [IS THIS ONE NEEDED GIVEN THE OTHERS?]

        # AEROSOLS
        # SW
        'tauaer_sw': [climt_inputs['tauaer_sw']], # Aerosol optical depth (iaer=10 only), Dimensions: (ncol,nlay,nbndsw), (non-delta scaled)
        'ssaaer_sw': [climt_inputs['ssaaer_sw']], # Aerosol single scattering albedo (iaer=10 only), Dimensions: (ncol,nlay,nbndsw), (non-delta scaled)      
        'asmaer_sw': [climt_inputs['asmaer_sw']], # Aerosol asymmetry parameter (iaer=10 only), Dimensions: (ncol,nlay,nbndsw), (non-delta scaled)
        'ecaer_sw': [None], # Aerosol optical depth at 0.55 micron (iaer=6 only), Dimensions: (ncol,nlay,naerec), (non-delta scaled)
        'tauaer_lw': [climt_inputs['tauaer_lw']]
    }
    import pdb; pdb.set_trace()
    return _rrtm_radiation_fortran.driver(*args)
