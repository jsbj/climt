#!/usr/bin/python

# RRTM Reference:
# Mlawer, E.J., S.J. Taubman, P.D. Brown,  M.J. Iacono and 
# S.A. Clough: RRTM, a validated correlated-k model for the 
# longwave. J. Geophys. Res., 102, 16,663-16,682, 1997         

# Gotta have this to let the browser know it's json
print "Content-type: application/json\n\n";

import json, sys, re, os, climt
from numpy import e, linspace, log
from subprocess import call
from os import rename, chdir
from math import copysign, floor, log10
import numpy

json_input = {} # json.load(sys.stdin)
if 'asdir' in json_input:
    sys.stderr.write(str(json_input['asdir']))
    json_input['asdif'] = json_input['asdir']
    json_input['aldir'] = json_input['asdir']
    json_input['aldif'] = json_input['asdir']
    

defaults = json.loads('{"ps": 1013, "scon": 1371, "asdir": 0.3, "asdif": 0.3, "aldir": 0.3, "aldif": 0.3, "tauaer_sw": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "ssaaer_sw": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "asmaer_sw": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "tauaer_lw": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "r_liq": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "r_ice": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "clwp": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "ciwp": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "cldf": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "ccl4": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "cfc11": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "cfc12": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "cfc22": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "co2": [355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355, 355], "n2o": [320.14772999999997, 320.14808, 320.12952, 320.17348, 320.20259, 320.2492, 320.18053000000003, 320.15103, 320.06952, 319.64703000000003, 317.94278, 314.85408, 309.93951, 302.8905, 297.28139, 293.30249, 286.54488999999995, 279.02988, 269.73828, 254.67132999999998, 231.32466, 199.50788999999997, 169.08091, 139.91885000000002, 117.2268, 103.31899, 94.382699, 87.561951, 82.404142, 75.596006, 66.9516, 54.150636000000006, 42.426844, 32.571123, 24.015852, 17.783966000000003, 12.921510000000001, 9.3075085, 6.6677854, 3.5912390999999997, 2.0309472, 1.7047587, 1.4732259, 1.3152129, 1.2046001, 1.1028871, 1.0173566, 0.9552473300000001, 0.9000983300000001, 0.8477577, 0.8001817499999999], "altitude": [1.1, 2.1, 2.9, 3.7, 4.5, 5.3, 6.1, 6.8, 7.5, 8.2, 8.9, 9.6, 10.3, 11, 11.7, 12.4, 13.6, 14.7, 15.8, 16.9, 18, 19.1, 20.2, 21.3, 22.5, 23.7, 24.9, 26.2, 27.6, 29.1, 30.8, 32.9, 34.9, 36.9, 38.9, 40.9, 42.9, 44.9, 47.2, 53.9, 56.4, 58.4, 60.3, 61.7, 63.1, 64.5, 65.7, 66.8, 67.9, 69, 70], "zen": 62, "Ts": 294.2, "h2o": [0.015946558, 0.011230157, 0.0076751928, 0.0052688639, 0.0036297729, 0.0023900282, 0.0017066754, 0.0012718296, 0.00095655693, 0.00069666741, 0.00050829613, 0.00036584702, 0.00024977655, 0.00013636267, 6.5472166e-05, 2.8419665e-05, 9.6973117e-06, 4.8207025e-06, 3.4318521e-06, 3.2663258e-06, 3.178493e-06, 3.1768304e-06, 3.2639416e-06, 3.4095149e-06, 3.5909502e-06, 3.8500998e-06, 4.0575464e-06, 4.251363e-06, 4.3863338e-06, 4.5309193e-06, 4.6839027e-06, 4.806785e-06, 4.9039072e-06, 4.9670398e-06, 5.016137e-06, 5.1013058e-06, 5.2471341e-06, 5.3810127e-06, 5.4697343e-06, 5.4735615e-06, 5.332653e-06, 5.1831207e-06, 5.0460312e-06, 4.8780507e-06, 4.7075605e-06, 4.5413699e-06, 4.3837813e-06, 4.2189254e-06, 4.0623413e-06, 3.9098322e-06, 3.7676771e-06], "p": [952.1147, 841.897, 755.3917, 685.0609, 620.7571, 561.5159, 506.7787, 458.9778, 417.8179, 379.9846, 345.1331, 313, 283.2681, 255.9648, 230.793, 207.5901, 179.6777, 149.8259, 125.467, 105.5072, 88.85838, 74.81903, 63.06029, 53.19867, 44.59128, 37.16316, 30.91292, 25.6397, 20.97451, 16.9346, 13.41941, 10.30125, 7.703475, 5.824757, 4.442682, 3.407392, 2.627624, 2.037819, 1.56118, 0.9634139, 0.5106084, 0.3820259, 0.2975729, 0.2388066, 0.1978831, 0.1639725, 0.1372726, 0.1161604, 0.098930135, 0.084255733, 0.072246842], "Tbound": [289.25, 284.6, 279.8, 275, 270.2, 265.4, 260.55, 256, 251.45, 246.9, 242.35, 237.86, 233.35, 228.8, 224.25, 219.7, 215.74, 215.7, 215.7, 215.7, 216.8, 218.03, 219.44, 220.76, 222.2, 223.57, 224.98, 226.71, 228.66, 231.81, 235.4, 239.99, 244.95, 249.84, 254.77, 259.73, 264.69, 269.65, 274.56, 270.71, 265.88, 261, 256.08, 251.32, 246.56, 241.8, 237.02, 232.18, 227.34, 222.5, 218.1], "ch4": [1700.7853, 1700.7861, 1700.6882, 1700.0174, 1696.7191, 1689.0904999999998, 1677.4702, 1662.5031999999999, 1646.9684, 1632.9801, 1622.3284999999998, 1607.1415, 1582.0669, 1556.2247, 1531.3253, 1508.0506, 1480.6419, 1447.9623, 1415.2675000000002, 1379.503, 1342.601, 1301.4651999999999, 1245.1943, 1172.2138, 1075.8682999999999, 965.1576, 854.01462, 771.0717099999999, 725.38978, 680.32085, 634.01592, 579.41355, 527.36578, 481.60666, 437.54815, 394.57359, 352.15132, 310.31249, 267.31394, 200.8872, 158.78383, 154.0019, 151.14806000000002, 150.15239, 150.18485, 150.16241, 150.13467, 150.23033, 150.28188, 150.26681, 150.18884], "T": [291.77, 287.03, 282.23, 277.43, 272.63, 267.83, 263.03, 258.3, 253.75, 249.2, 244.65, 240.13, 235.64, 231.1, 226.55, 222.01, 216.81, 215.71, 215.7, 215.7, 216.18, 217.39, 218.72, 220.08, 221.46, 222.88, 224.24, 225.81, 227.61, 230.17, 233.52, 237.51, 242.34, 247.27, 252.17, 257.13, 262.09, 267.05, 272, 274.41, 268.77, 263.53, 258.75, 253.76, 249, 244.24, 239.61, 234.65, 229.81, 224.97, 220.34], "co": [1.4735235e-07, 1.4203219e-07, 1.3746356e-07, 1.338817e-07, 1.3135738e-07, 1.3046302e-07, 1.293139e-07, 1.2701938e-07, 1.2377659e-07, 1.1940332e-07, 1.1352941e-07, 1.0700342e-07, 1.0015444e-07, 9.3152551e-08, 8.5588468e-08, 7.7191764e-08, 6.3881643e-08, 4.8797485e-08, 3.7298612e-08, 2.8723687e-08, 2.2545748e-08, 1.7379815e-08, 1.4111547e-08, 1.2622904e-08, 1.2397807e-08, 1.3167179e-08, 1.4350868e-08, 1.5625453e-08, 1.6708778e-08, 1.8091109e-08, 1.9843396e-08, 2.1874927e-08, 2.384691e-08, 2.5646894e-08, 2.7513584e-08, 2.9431952e-08, 3.0938047e-08, 3.230932e-08, 3.3800561e-08, 3.6464382e-08, 3.9601694e-08, 4.2654523e-08, 4.5695458e-08, 4.9774858e-08, 5.4377978e-08, 5.9385144e-08, 6.5223382e-08, 7.4618846e-08, 8.5339593e-08, 9.7556516e-08, 1.1081534e-07], "o3": [0.031872162, 0.035456235, 0.039477314, 0.043921091, 0.048850309999999994, 0.054422609999999996, 0.061250461, 0.069855773, 0.079463597, 0.08915115, 0.10168034, 0.1155858, 0.13068458, 0.16048106, 0.19350828, 0.22751290999999998, 0.304286, 0.43981947, 0.52382995, 0.63216254, 0.82302279, 1.2512421999999999, 1.8039109, 2.2908109, 2.8324889, 3.4517834, 4.2219771999999995, 5.032683899999999, 5.6775239, 6.3139009, 6.9619100000000005, 7.772886399999999, 8.524654700000001, 8.8305105, 8.490472299999999, 7.5621829, 6.2966351, 5.1043844, 4.0821087, 2.8155102000000003, 1.803627, 1.545081, 1.3594723, 1.1832445999999999, 1.0330702, 0.90162695, 0.78788491, 0.67509507, 0.57978644, 0.49771251, 0.42984522000000003], "o2": [0.20897518, 0.20897572, 0.2089678, 0.2089866, 0.20899189, 0.20899543, 0.20899996, 0.20900373, 0.20900458, 0.20900519, 0.20900649, 0.20900634, 0.20900698, 0.20900562, 0.20900711, 0.20900925, 0.20900522, 0.20899965, 0.20899954, 0.20899963, 0.20899959, 0.20899966, 0.20899986, 0.20899987, 0.20900002, 0.20899989, 0.20899986, 0.2090022, 0.20900251, 0.2090067, 0.2090057, 0.20900536, 0.20900574, 0.20900482, 0.20900646, 0.20900702, 0.20900613, 0.20900463, 0.2090015, 0.20900197, 0.20901358, 0.2090466, 0.20902328, 0.20906644, 0.20911193, 0.20908101, 0.20904104, 0.20916539, 0.20922786, 0.20919746, 0.20908001], "lev": [891.46, 792.287, 718.704, 651.552, 589.841, 532.986, 480.526, 437.556, 398.085, 361.862, 328.507, 297.469, 269.015, 243, 218.668, 196.44, 162.913, 136.511, 114.564, 96.4903, 81.2, 68.4286, 57.6936, 48.6904, 40.5354, 33.733, 28.1201, 23.1557, 18.7914, 15.0693, 11.8006, 8.78628, 6.61328, 5.03469, 3.85333, 2.96408, 2.2918, 1.78227, 1.339, 0.589399, 0.430705, 0.333645, 0.261262, 0.216491, 0.179393, 0.148652, 0.1255, 0.106885, 0.091031, 0.077529, 0.067]}')

# json_input = dict(json_input.items() + defaults.items())
new_json_input = {}
for key in defaults:
    new_json_input[key] = defaults[key]
for key in json_input:
    new_json_input[key] = json_input[key]

json_input = new_json_input

# import sys; sys.stderr.write(str(json_input['cfc11']))
if 'co2' in json_input:
    json_input['co2'] = [co2 / 1.e6 for co2 in json_input['co2']]
if 'ch4' in json_input:
    json_input['ch4'] = [ch4 / 1.e9 for ch4 in json_input['ch4']]
if 'n2o' in json_input:
    json_input['n2o'] = [n2o / 1.e9 for n2o in json_input['n2o']]
if 'o3' in json_input:
    json_input['o3'] = [o3 / 1.e6 for o3 in json_input['o3']]
if 'cfc11' in json_input:
    json_input['cfc11'] = [cfc11 / 1.e12 for cfc11 in json_input['cfc11']]
    json_input['cfc12'] = [cfc12 / 1.e12 for cfc12 in json_input['cfc12']]
    json_input['cfc22'] = [cfc22 / 1.e12 for cfc22 in json_input['cfc22']]    
    json_input['ccl4'] = [ccl4 / 1.e12 for ccl4 in json_input['ccl4']]
r = climt.radiation(scheme='rrtm', **json_input)

# In case it's needed:
def sigdig(x, digits=1):
    if x:
        x = copysign(round(x, -int(floor(log10(abs(x)))) + (digits - 1)), x)
    return x

results = {
    'swuflx': [sigdig(float(f), 3) for f in r['swuflx']],
    'swdflx': [sigdig(float(f), 3) for f in r['swdflx']],
    'lwuflx': [sigdig(float(f), 3) for f in r['lwuflx']],
    'lwdflx': [sigdig(float(f), 3) for f in r['lwdflx']],
    'uflx': [sigdig(float(r['swuflx'][i] + r['lwuflx'][i]), 3) for i in range(len(r['swuflx']))],
    'dflx': [sigdig(float(r['swdflx'][i] + r['lwdflx'][i]), 3) for i in range(len(r['swdflx']))],
    'LwToa': sigdig(float(r['LwToa']), 3),
    'SwToa': sigdig(float(r['SwToa']), 3),
    'net_toa': round(sigdig(float(r['LwToa']) + float(r['SwToa']), 3))
}
output = {}

for key in new_json_input:
    output[key] = new_json_input[key]
for key in results:
    output[key] = results[key]
output['co2'] = [co2 * 1.e6 for co2 in output['co2']]
output['ch4'] = [ch4 * 1.e9 for ch4 in output['ch4']]
output['n2o'] = [n2o * 1.e9 for n2o in output['n2o']]
output['o3'] = [o3 * 1.e6 for o3 in output['o3']]
output['cfc11'] = [cfc11 * 1.e12 for cfc11 in output['cfc11']]
output['cfc12'] = [cfc12 * 1.e12 for cfc12 in output['cfc12']]
output['cfc22'] = [cfc22 * 1.e12 for cfc22 in output['cfc22']]
output['ccl4'] = [ccl4 * 1.e12 for ccl4 in output['ccl4']]
print json.dumps(output)


def load_atmosphere(atmosphere):
    # loads a dictionary from a data file stored in /atmospheres
    f = open('atmospheres/' + atmosphere + '.json', 'r')
    atmosphere = json.load(f)
    f.close()
    
    return atmosphere