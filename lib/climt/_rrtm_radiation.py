from scipy.interpolate import interp1d
from math import cos, pi
import _rrtm_radiation_fortran
from numpy import ndarray
INPUTS = [
            # 'do_sw', #     0  Shortwave switch     (integer)       1   1 / 0 => do / do not compute SW        
            # 'do_lw', #     0  Longwave switch      (integer)       1   1 / 0 => do / do not compute LW       
            'p', #       1-3  Atmospheric pressure     mb               Default is equispaced 0-ps. p[0] is top level
            'lev',
            'T', #       1-3  Temperature              K        283.15  Isothermal
            'Tbound',
            'ps', #      0-2  Surface pressure         mb       1000.
            'Ts', #      0-2  Surface temperature      K        283.15
            'q', #       1-3  Specific humidity        g/kg     1.e-5
            'h2o',
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
            'tauaer_lw', # Aerosol optical depth (iaer=10 only), Dimensions: (ncol,nlay,nbndlw), (non-delta scaled)
            'Cpd',
            'tauc_lw'
]

def driver(*args):
    # wavenumber bands used by RRTM:
    SW_BANDS = range(14)
    LW_BANDS = range(16)
    
    # gotta translate between the APIs:
    climt_inputs = dict(zip(INPUTS, args))

    number_of_layers = len(climt_inputs['T'])
    
    climt_inputs['pbound'] = climt_inputs['ps'][0].tolist() + climt_inputs['lev'].tolist()
    climt_inputs['Tbound'] = [climt_inputs['Ts'][0][0]] + [T[0][0] for T in climt_inputs['Tbound']]
    interpolated_p = interp1d(range(number_of_layers + 1), climt_inputs['pbound'])
    interpolated_T = interp1d(range(number_of_layers + 1), climt_inputs['Tbound'])
    
    clouds = 1 if 'cldf' in climt_inputs else 0
    
    # import sys; sys.stderr.write(str(climt_inputs['cldf']))
    
    for key in ['co2', 'ch4', 'n2o', 'o2', 'cfc11', 'cfc12', 'cfc22', 'ccl4']:
        if not hasattr(climt_inputs[key], '__iter__'):
            climt_inputs[key] = [climt_inputs[key]] * number_of_layers
    for key in ['lw_surface_emissivity']:
        if not hasattr(climt_inputs[key], '__iter__'):
            climt_inputs[key] = [climt_inputs[key]] * len(LW_BANDS)
    
    if 'h2o' in climt_inputs:
        h2o_concentration = [[h2o[0][0] for h2o in climt_inputs['h2o']]]
    else:
        h2o_concentration = [[(((q/1000.)/(1. - (q/1000.)))*1.607793)[0][0] for q in climt_inputs['q']]]
    
    for key in ['tauaer_sw', 'ssaaer_sw', 'asmaer_sw', 'tauaer_lw']:
        if not hasattr(climt_inputs[key][0], '__iter__'):
            climt_inputs[key] = [[value] * len(locals()[key[-2:].upper() + '_BANDS']) for value in climt_inputs[key]]
    

    
    rrtm_inputs = [
        # len(LW_BANDS), # 'nbndlw' - ideally, these four variables would be called with us,
        140, # 'ngptlw',
        # len(SW_BANDS), # 'nbndsw'
        112, # 'ngptsw'
        # GENERAL, used in both SW and LW
        1, #'iplon':  # index of the column (1 indexed for Fortran)
        # 1, #'ncol':  # number of columns
        # number_of_layers, #'nlay': 
        clouds, # Cloud overlap method, 0: Clear only, 1: Random, 2: Maximum/random, #'icld':  3: Maximum
        150, #'permuteseed_sw':  # used for monte carlo clouds; must differ from permuteseed_lw by number of subcolumns
        300, #'permuteseed_lw':  # learn about these later...
        1, #'irng': 
        0, #'idrv':  # whether to also calculate the derivative of flux with respect to surface temp
        climt_inputs['Cpd'], #'cpdair':  
        [[interpolated_p(i + .5) for i in range(number_of_layers)]], #'play':  # pressure in each layer
        [climt_inputs['pbound']], # pressure at boundaries of each layer, #'plev':  with a linear extrapolation for the "top"
        [[interpolated_T(i + .5) for i in range(number_of_layers)]], #'tlay':  # temperature in each layer
        [climt_inputs['Tbound']],
        [climt_inputs['Ts']], #'tsfc': 
        # GASES, used in both SW and LW
        h2o_concentration, #'h2ovmr':  # convert from g/kg to volume mixing ration using molecular weight of dry air / water vapor
        [[o3[0][0] * 0.603428 for o3 in climt_inputs['o3']]], #'o3vmr':  # convert from kg/kg to volume mixing ratio using molecular weight of dry air / ozone
        [[co2 for co2 in climt_inputs['co2']]], #'co2vmr': 
        [[ch4 for ch4 in climt_inputs['ch4']]], #'ch4vmr': 
        [[n2o for n2o in climt_inputs['n2o']]], #'n2ovmr': 
        [climt_inputs['o2']], #'o2vmr': 
        [[0. or cfc11 for cfc11 in climt_inputs['cfc11']]], #'cfc11vmr': 
        [[0. or cfc12 for cfc12 in climt_inputs['cfc12']]], #'cfc12vmr': 
        [[0. or cfc22 for cfc22 in climt_inputs['cfc22']]], #'cfc22vmr': 
        [climt_inputs['ccl4']], #'ccl4vmr': 
        # SURFACE OPTICAL PROPERTIES
        # SW
        [climt_inputs['aldif'][0][0]], #'aldif': 
        [climt_inputs['aldir'][0][0]], #'aldir': 
        [climt_inputs['asdif'][0][0]], #'asdif': 
        [climt_inputs['asdir'][0][0]], #'asdir': 
        # LW
        [[1. or 1 - emis for emis in climt_inputs['lw_surface_emissivity']]], #'emis': 
        # THE SUN - SW
        [cos(climt_inputs['zen'][0][0] * 2 * pi / 360.)], #'coszen':  # cosine of the solar zenith angle
        1., #'adjes':  # flux adjustment for earth/sun distance (if not dyofyr)
        0, # day of the year, #'dyofyr':  used to get Earth/Sun distance (if not adjes)
        climt_inputs['scon'], #'scon':  # solar constant
        # CLOUDS, SW see http://www.arm.gov/publications/proceedings/conf16/extended_abs/iacono_mj.pdf
        2, #'inflgsw':  # Flag for cloud optical properties
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
        2, #'inflglw': 
        1, #'iceflgsw':  # Flag for ice particle specification
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
        1, #'iceflgslw': 
        1, #'liqflgsw':  # Flag for liquid droplet specification
            # LIQFLAG = 0 the optical depths (gray) due to water clouds are computed as in CCM3.
            #         = 1 the optical depths (non-gray) due to water clouds are computed by a method
            #             based on the parameterization of water clouds due to Y.X. Hu and K. Stamnes,
            #             J. Clim., 6, 728-742 (1993).
        1, #'liqflglw': 
        [[[0.]* number_of_layers]] * len(SW_BANDS), #'tauc_sw':  # In-cloud optical depth [IS THIS ONE NEEDED GIVEN THE OTHERS?]
        [[climt_inputs['tauc_lw'] or [0.] * number_of_layers]] * len(LW_BANDS), #'tauc_lw':  # in-cloud optical depth 
        [[c[0][0] for c in climt_inputs['cldf']]], #'cldfrac':  # layer cloud fraction
        [[climt_inputs['cloud_single_scattering_albedo'] or [0.] * number_of_layers]] * len(SW_BANDS), #'ssac_sw':  # In-cloud single scattering albedo
        [[climt_inputs['cloud_asymmetry_parameter'] or [0.] * number_of_layers]] * len(SW_BANDS), #'asmc_sw':  # In-cloud asymmetry parameter
        [[climt_inputs['cloud_forward_scattering_fraction'] or [0.] * number_of_layers]] * len(SW_BANDS), #'fsfc_sw':  # In-cloud forward scattering fraction (delta function pointing forward "forward peaked scattering")
        [[c[0][0] for c in climt_inputs['ciwp']]], #'ciwp':  # in-cloud ice water path (g/m2)
        [[c[0][0] for c in climt_inputs['clwp']]], #'clwp':  # in-cloud liquid water path (g/m2)
        [[c[0][0] for c in climt_inputs['r_ice']]], #'reic':  # Cloud ice particle effective size (microns)
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
        [[c[0][0] for c in climt_inputs['r_liq']]], #'relq':  # Cloud water drop effective radius (microns)


        # AEROSOLS
        # SW
        [climt_inputs['tauaer_sw'] or [[0.] * len(SW_BANDS)] * number_of_layers], # Aerosol optical depth (iaer=10 only), Dimensions: (ncol,nlay,nbndsw), #'tauaer_sw':  (non-delta scaled)
        [climt_inputs['ssaaer_sw'] or [[0.] * len(SW_BANDS)] * number_of_layers], # Aerosol single scattering albedo (iaer=10 only), Dimensions: (ncol,nlay,nbndsw), #'ssaaer_sw':  (non-delta scaled)      
        [climt_inputs['asmaer_sw'] or [[0.] * len(SW_BANDS)] * number_of_layers], # Aerosol asymmetry parameter (iaer=10 only), Dimensions: (ncol,nlay,nbndsw), #'asmaer_sw':  (non-delta scaled)
        [[[0.] * 6] * number_of_layers], # Aerosol optical depth at 0.55 micron (iaer=6 only), Dimensions: (ncol,nlay,naerec), #'ecaer_sw':  (non-delta scaled)
        [climt_inputs['tauaer_lw'] or [[0.] * len(LW_BANDS)] * number_of_layers] #'tauaer_lw': 
    ]
    # import pdb; pdb.set_trace()
    output = _rrtm_radiation_fortran.driver(*rrtm_inputs)

    # 
    new_output = (
        output[0][0], # swuflx
        output[1][0], # swdflx
        output[6][0], # lwuflx
        output[7][0],  # lwdflx
        output[1][0][-1] - output[0][0][-1], # swToA
        output[7][0][-1] - output[6][0][-1] # lwToA
    )
    return new_output