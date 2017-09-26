"""
@file irfrun.py

@brief Run IRF generation (makeirf,makefits,livetime_corr) in the
target directory.
"""

import sys
import os
import glob
import re
import argparse
import yaml

if not 'INST_DIR' in os.environ:
    raise Exception('ST installation directory not defined.')
else:
    inst_dir = os.environ['INST_DIR']

sys.path.append(os.path.join(inst_dir,'irfs/handoff_response/python'))

from run_livetime_corr import corr_fit
from irfutils import get_class_types_from_xml

def makeirf(irfname, ext='', overwrite=False):

    setup_file = 'setup%s.py'%ext
    skim_file = 'skim_%s%s.root'%(irfname,ext)
    param_file = 'par_%s%s.root'%(irfname,ext)

    if not os.path.isfile(setup_file):
        raise Exception('No setup.py in target directory.')
    
    #if not os.path.isfile(skim_file) or overwrite:
    #    cmd = 'prune %s'%setup_file
    #    print cmd
    #    os.system(cmd)

    if not os.path.isfile(param_file) or overwrite:
        cmd = 'makeirf %s'%setup_file
        print cmd
        os.system(cmd)

    cmd = 'makefits %s %s'%(param_file,irfname)
    print cmd
    os.system(cmd)

usage = "%(prog)s [options] [target_dir]"
description = """This utility runs IRF generation in the target directory."""
parser = argparse.ArgumentParser(usage=usage, description=description)

parser.add_argument('target_dir', nargs='*',
                    help='IRF target directories.')

parser.add_argument('--classdefs',
                    default=None, 
                    help='Set the XML class definitions file.  This will be used to determine '
                    'the set of classes and event types for which IRFs will be generated.')

parser.add_argument('--batch',default=False,action='store_true')

parser.add_argument('--overwrite',default=False,action='store_true')

parser.add_argument('--no_livetime',default=False,action='store_true',
                    help='Skip generation of livetime correction.')

#parser.add_argument('--config',
#                    default = None,required=True,
#                    help = 'Set the configuration file for IRF generation.  '
#                    'This file defines scaling and fit parameters as well as '
#                    'the energy range and number of thrown events in each '
#                    'MC sample.')

args = parser.parse_args()

if args.classdefs is not None:

    class_version, classnames, typenames = get_class_types_from_xml(
        args.classdefs)
    for classname in classnames:
        for typename in typenames:
            irfname = '%s_%s_%s' % (class_version, classname, typename)
            args.target_dir += [irfname]

args.target_dir = [os.path.abspath(t) for t in args.target_dir]
cwd = os.getcwd()

if args.batch:

    for d in args.target_dir:

        if not os.path.isdir(d):
            continue
        
        logfile = os.path.join(d,'irfrun.log')
        shfile = os.path.join(d,'irfrun.sh')
        exe = os.path.join(os.path.dirname(__file__),'irfrun.py')
        bash_script = "python {exe} %s"%(d)
        if args.no_livetime:
            bash_script += ' --no_livetime '
        if args.overwrite:
            bash_script += ' --overwrite '            
        scriptexe = 'irfrun.sh'
        with open(shfile, 'wt') as f:
            f.write(bash_script.format(exe=exe))

        cmd = 'bsub -oo %s -W 500 bash %s'%(logfile, shfile)
        print(cmd)
        os.system(cmd)
        
    sys.exit(0)
    
for target_dir in args.target_dir:

    if not os.path.isfile(os.path.join(target_dir,'setup.py')):
        raise Exception('No setup.py in target directory.')

    if not os.path.isfile(os.path.join(target_dir,'config.yaml')):
        raise Exception('No config.yaml in target directory.')
    
    #config = yaml.load(open(args.config,'r'))
    config = yaml.load(open(os.path.join(target_dir,'config.yaml'),'r'))

    os.chdir(target_dir)

    # Generate IRFs     
    makeirf(config['irf_name'],overwrite=args.overwrite)

    if args.no_livetime or not 'Livetime' in config:
        continue

    # Generate IRFs for each livetime bin
    livetime_bins = sorted(config['Livetime']['MC'].keys())
    for ltb in livetime_bins:
        makeirf(config['irf_name'],'_%s'%(ltb),
                overwrite=args.overwrite)

    # Generate the livetime correction
    lts = config['Livetime']['lts']
    avltf=config['Livetime']['avltf']
    path='./'
    reffile=path+'aeff_'+config['irf_name']+'.fits'

    cf=corr_fit(reffile,avltf,config['irf_name'],
                min_ebin=8,max_ebin=68)

    for lt, ltb in zip(lts,livetime_bins):
        fn=("aeff_%s_%s.fits") % (config['irf_name'],ltb)
        cf.add_aeff_file([fn,lt],[25.e6,1.25,5.75])

    old_p0=config['Livetime']['p0_pars']

    #do the fit
    cf.fit(old_p0,fits_output="")
    cf.writein(reffile)


os.chdir(cwd)
