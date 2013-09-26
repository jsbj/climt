#!/usr/bin/env python

import os,glob,string,sys
from numpy.distutils.core import setup, Extension
from numpy.distutils import fcompiler
from distutils.dep_util import newer

## -------- set these
KM = 26
JM = 1
IM = 1
NC_INC = '/usr/local/include'
NC_LIB = '/usr/local/lib'
##----------------------

if '--lite' in sys.argv:
    sys.argv.pop(sys.argv.index('--lite'))
    Lite = True
else:
    Lite = False

Extensions = [
    {'name':'grid',
     'dir':'src/grid'},
    {'name':'timestep',
     'dir':'src/timestep'},
#    {'name':'thermodyn',
#     'dir':'src/thermodyn'},
#    {'name':'emanuel_convection',
#     'dir':'src/convection/emanuel'},
#    {'name':'hard_adjustment',
#     'dir':'src/convection/hard'},
#    {'name':'sbm_convection',
#     'dir':'src/convection/hard'},
#    {'name':'axisymmetric_dynamics',
#     'dir':'src/dynamics/axisymmetric'},
#    {'name':'slab_ocean',
#     'dir':'src/ocean/slab_ocean'},
#    {'name':'ccm3_radiation',
#     'dir':'src/radiation/ccm3',
#     'cppflags':'-DSUN -DPLON=%i -DPLEV=%i -DPLEVR=%i' % (IM,KM,KM)},
#    {'name':'cam3_radiation',
#     'dir':'src/radiation/cam3',
#     'cppflags':'-DPLEV=%i' % KM,
#     'lib':['netcdf','netcdff'],
#     'libdir': [NC_LIB],
#     'incdir': [NC_INC]},
#    {'name':'chou_radiation',
#     'dir':'src/radiation/chou'},
#    {'name':'greygas_radiation',
#     'dir':'src/radiation/greygas'},
#    {'name':'ozone',
#     'dir':'src/radiation/ozone'},
    {'name':'insolation',
     'dir':'src/radiation/insolation'},
#    {'name':'ccm3_turbulence',
#     'dir':'src/turbulence/ccm3',
#     'cppflags':'-DPLON=%i -DPLEV=%i' % (IM,KM)},
#    {'name':'simple_turbulence',
#     'dir':'src/turbulence/simple'},
    {'name':'rrtm_radiation_fortran',
     'dir':'src/radiation/rrtm'}
    ]

# define extensions that will be built when the --lite option is used
LiteExtensionsNames = ['grid','timestep','insolation','ozone','thermodyn','ccm3_radiation']
ExtensionsLite = []
for ext in Extensions:
    if ext['name'] in LiteExtensionsNames:
        ExtensionsLite.append(ext)

# figure out which compiler we're goint to use
compiler = fcompiler.get_default_fcompiler()
for i in range(len(sys.argv)):
    if '--fcompiler' in sys.argv[i]:
        compiler = sys.argv.pop(i)
        compiler = compiler[compiler.index('=')+1:]
print 'Using %s compiler' % compiler

# set some fortran compiler-dependent flags
if compiler == 'gnu95':
    compiler = 'gfortran'
    f77flags='-ffixed-line-length-132 -fdefault-real-8'
    f90flags='-fdefault-real-8'
elif compiler == 'intel' or compiler == 'intelem':
    f77flags='-132 -r8 -w95 -w90 -mp'
    f90flags='-r8 -w95 -mp'
elif compiler == 'ibm':
    f77flags='-qautodbl=dbl4 -qsuffix=f=f:cpp=F -qfixed=132'
    f90flags='-qautodbl=dbl4 -qsuffix=f=f90:cpp=F90 -qfree=f90'
else:
    print 'Sorry, compiler %s not supported' % compiler

for ExtList in [Extensions,ExtensionsLite]:
    for i in range(len(ExtList)):
        Defaults = {'cppflags':'-DIM=%i -DJM=%i -DKM=%i' % (IM,JM,KM),
                    'f77flags':f77flags,
                    'f90flags':f90flags}
        Defaults.update(ExtList[i])
        ExtList[i] = Defaults
    if compiler == 'ibm':
        for ext in ExtList:
            ext['cppflags']='-WF,'+string.join(ext['cppflags'].split(),',')

def getSources(dir):
    #Gets list of source files for extensions
    SrcFile = os.path.join(dir,'sources_in_order_of_compilation')
    if os.path.exists(SrcFile):
        Sources = open(SrcFile).readlines()
        Sources = [os.path.join(dir,s[:-1]) for s in Sources]
    else:
        Sources = []
        for pattern in ['*.f','*.F','*.f90','*.F90']:
            Sources += glob.glob(os.path.join(dir,pattern))
            Sources += glob.glob(os.path.join(dir,'src',pattern))        
    return Sources

def buildNeeded(target,src):
    #Checks if source code is newer than extension, so extension needs to be rebuilt
    target = os.path.join('lib/climt',target)
    if not os.path.exists(target):
        return True
    for file in src:
        if newer(file,target):
            return True
    print 'Extension %s is up to date' % os.path.basename(target)
    return False

def build_ext(name=None, dir=None, cppflags='', f77flags='', f90flags='', \
              lib='', libdir='', incdir=''):
    if name == 'rrtm_radiation_fortran':
	    #Builds an extension
	    os.system('rm -rf tmp; mkdir tmp')
	    src = getSources(dir)
	    target = '_%s.so' % name
	    driver = glob.glob(os.path.join(dir,'Driver.f*'))[0]
	    f77flags = '-c %s %s' % (cppflags,f77flags)
	    f90flags = '-c -fPIC -fno-range-check %s %s' % (cppflags,f90flags)
	    if buildNeeded(target,src):
		print '\n Building %s ... \n' % os.path.basename(target)
		for filename in src:
		    os.system('cp %s tmp/' % filename)
		f77src = [filename.split('/')[-1] for filename in src if filename[-2:] in ['.f', '.F']]
		f90src = [filename.split('/')[-1] for filename in src if filename[-4:] in ['.f90', '.F90']]
		if len(f77src) > 0:
		    f77command = 'cd tmp;' + ' '.join((compiler, f77flags, ' '.join(f77src)))
		    print f77command
		    os.system(f77command)
		    
		if len(f90src) > 0:
		    f90command = 'cd tmp;' + ' '.join((compiler, f90flags, ' '.join(f90src)))
		    print f90command
		    os.system(f90command)
		os.system('cp .f2py_f2cmap tmp/')    
		os.system('cd tmp; f2py -m _%s -h _%s.pyf --overwrite-signature %s' % (name,name,driver.split('/')[-1]))
		os.system('cd tmp; f2py -c --fcompiler=gnu95 _%s.pyf --build-dir .  *.o' % (name))
		os.system('cd tmp; gcc -pthread -shared ./src.linux-x86_64-2.6/_rrtm_radiation_fortranmodule.o ./src.linux-x86_64-2.6/fortranobject.o Driver.o mcica_random_numbers.o mcica_subcol_gen_lw.o mcica_subcol_gen_sw.o parkind.o parrrsw.o parrrtm.o rrlw_cld.o rrlw_con.o rrlw_kg01.o rrlw_kg02.o rrlw_kg03.o rrlw_kg04.o rrlw_kg05.o rrlw_kg06.o rrlw_kg07.o rrlw_kg08.o rrlw_kg09.o rrlw_kg10.o rrlw_kg11.o rrlw_kg12.o rrlw_kg13.o rrlw_kg14.o rrlw_kg15.o rrlw_kg16.o rrlw_ncpar.o rrlw_ref.o rrlw_tbl.o rrlw_vsn.o rrlw_wvn.o rrsw_aer.o rrsw_cld.o rrsw_con.o rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o rrsw_kg28.o rrsw_kg29.o rrsw_ncpar.o rrsw_ref.o rrsw_tbl.o rrsw_vsn.o rrsw_wvn.o rrtmg_lw_cldprmc.o rrtmg_lw_init.o rrtmg_lw_k_g.o rrtmg_lw_rad.o rrtmg_lw_rtrnmc.o rrtmg_lw_setcoef.o rrtmg_lw_taumol.o rrtmg_sw_cldprmc.o rrtmg_sw_init.o rrtmg_sw_k_g.o rrtmg_sw_rad.o rrtmg_sw_reftra.o rrtmg_sw_setcoef.o rrtmg_sw_spcvmc.o rrtmg_sw_taumol.o rrtmg_sw_vrtqdr.o -L/usr/lib64 -lpython2.6 -lgfortran -o ./_rrtm_radiation_fortran.so')
		os.system('mv tmp/_%s.so lib/climt' % name)
		os.system('rm -rf tmp')
		# # generate signature file
		# os.system('f2py --overwrite-signature %s -m _%s -h _%s.pyf'%(driver,name,name))
		# # compile extension
		# F2pyCommand = []
		# F2pyCommand.append('f2py -c -m _%s' % name)
		# F2pyCommand.append('--fcompiler=%s' % compiler)
		# F2pyCommand.append('-I%s' % dir)
		# F2pyCommand.append('-I%s' % os.path.join(dir,'include'))
		# F2pyCommand.append('-I%s' % os.path.join(dir,'src'))
		# F2pyCommand.append('-I%s' % os.path.join(dir,'src','include'))
		# if incdir is not '':
		#     for i in incdir:
		#         F2pyCommand.append('-I%s' % i)
		# if libdir is not '':
		#     for i in libdir:
		#         F2pyCommand.append('-L%s' % i)
		# if lib is not '':
		#     for i in lib:
		#         F2pyCommand.append('-l%s' % i)
		# F2pyCommand.append('--f77flags=%s' % f77flags)
		# F2pyCommand.append('--f90flags=%s' % f90flags)
		# F2pyCommand.append('_%s.pyf' % name)
		# F2pyCommand.append('%s' % string.join(src))
		# F2pyCommand = string.join(F2pyCommand)
		# print F2pyCommand
		# import pdb;pdb.set_trace()
		# if os.system(F2pyCommand) > 0:
		#     print '+++ Compilation failed'
		#     sys.exit()
		# os.system('mv -f _%s.so lib/climt' % name)
		# os.system('rm -f _%s.pyf' % name)
    else:
            src = getSources(dir)
            target = '_%s.so' % name
	    driver = glob.glob(os.path.join(dir,'Driver.f*'))[0]
	    f77flags = '"%s %s"' % (cppflags,f77flags)
	    f90flags = '"%s %s"' % (cppflags,f90flags)
	    if buildNeeded(target,src):
		print '\n Building %s ... \n' % os.path.basename(target)
		# generate signature file
		os.system('f2py --overwrite-signature %s -m _%s -h _%s.pyf'%(driver,name,name))
		# compile extension
		F2pyCommand = []
		F2pyCommand.append('f2py -c -m _%s' % name)
		F2pyCommand.append('--fcompiler=%s' % compiler)
		F2pyCommand.append('-I%s' % dir)
		F2pyCommand.append('-I%s' % os.path.join(dir,'include'))
		F2pyCommand.append('-I%s' % os.path.join(dir,'src'))
		F2pyCommand.append('-I%s' % os.path.join(dir,'src','include'))
		if incdir is not '':
		    for i in incdir:
			F2pyCommand.append('-I%s' % i)
		if libdir is not '':
		    for i in libdir:
			F2pyCommand.append('-L%s' % i)
		if lib is not '':
		    for i in lib:
			F2pyCommand.append('-l%s' % i)
		F2pyCommand.append('--f77flags=%s' % f77flags)
		F2pyCommand.append('--f90flags=%s' % f90flags)
		F2pyCommand.append('_%s.pyf' % name)
		F2pyCommand.append('%s' % string.join(src))
		F2pyCommand = string.join(F2pyCommand)
		print F2pyCommand
		if os.system(F2pyCommand) > 0:
		    print '+++ Compilation failed'
		    sys.exit()
		os.system('mv -f _%s.so lib/climt' % name)
		os.system('rm -f _%s.pyf' % name)

def setupClimt():
    # Build all extensions
    for ext in Extensions: build_ext(**ext)

    # Finish the setup
    # note: setup() cannot copy directories, and falls over
    # trying to copy the CVS directory in climt/lib/data
    # workaround: make data list which specifically excludes CVS
    os.chdir('lib/climt')
    DataFiles = []
    for File in glob.glob('data/*/*'):
        if 'CVS' not in File:
            DataFiles.append(File)
    print DataFiles
    os.chdir('../..')
    
    setup(name = "CliMT",
          version = open('Version').read()[:-1],
          description = "Climate modelling and diagnostics toolkit",
          author = "Rodrigo Caballero",
          author_email = "rodrigo@misu.su.se",
          url = "http://people.su.se/~rcaba/climt",
          packages = ['climt'],
          package_dir = {'':'lib'},
          package_data = {'climt':['*.so']+DataFiles})


def setupClimtLite():
    # Build all extensions
    for ext in ExtensionsLite: build_ext(**ext)
    os.system('mkdir -p lib/climt_lite')
    ClimtLiteFiles = ['__init__.py', '__version__.py', '_ccm3_radiation.so', '_grid.so',
     '_insolation.so', '_ozone.so', '_thermodyn.so', '_timestep.so', 'component.py',
     'grid.py', 'insolation.py', 'mathutil.py', 'ozone.py', 'parameters.py',
     'radiation.py', 'state.py', 'thermodyn.py', 'utils.py', 'io.py', 'plot.py']
    for file in ClimtLiteFiles:
        os.system('cp lib/climt/%s lib/climt_lite/%s' % (file,file))
    setup(name         = "CliMT-lite",
          version      = open('Version').read()[:-1],
          description  = "Climate modelling and diagnostics toolkit, lite version",
          author       = "Rodrigo Caballero",
          author_email = "rodrigo@misu.su.se",
          url          = "http://people.su.se/~rcaba/climt",
          packages    = ['climt_lite'],
          package_dir = {'climt_lite':'lib/climt_lite'},
          package_data = {'climt_lite':['*.so']})

if Lite:
    setupClimtLite()
else:
    setupClimt()
