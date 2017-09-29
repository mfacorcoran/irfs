"""
@file irfsetup.py

@brief Generate python setup files and directory structure for
makeirf/makefits.
"""
from __future__ import absolute_import, division, print_function

import yaml
import re
import collections
import os
import sys
import glob
import array
from xml.dom import minidom
import xml.etree.cElementTree as ElementTree
import argparse
import itertools
import numpy as np
import ROOT

from irfutils import get_branches, get_cuts_from_xml, extract_generated_events
from irfutils import get_class_types_from_xml, makedir

def write_config(config, name):

    out_str = ''
    for k, v in sorted(config.items()):
        if isinstance(v, str):
            out_str += r"{}.{} = '{}'".format(name, k, v) + '\n'
        else:
            out_str += r"{}.{} = {}".format(name, k, v) + '\n'
    return out_str


def make_generated_events_str(config):
    """Generate configuration string for the number of generated
    events in each MC sample."""

    ngen = [v['ngen'] for k, v in config.items()]
    logemin = [v['logemin'] for k, v in config.items()]
    logemax = [v['logemax'] for k, v in config.items()]
    out_str = ''
    out_str += "Data.generated = %s\n" % str(ngen)
    out_str += "Data.logemin = %s\n" % str(logemin)
    out_str += "Data.logemax = %s\n" % str(logemax)
    return out_str


def makesetup(class_version, classname, variant, irfname, config, sel, livetime_bin=None):

    if livetime_bin is None:
        generated = config['MC']
    else:
        generated = {livetime_bin : config['Livetime']['MC'][livetime_bin] }

    setup_string = """
from IRFdefault import *
from math import *\n
"""

    setup_string += "Prune.fileName = ''\n"    
    #if livetime_bin is None:
    #    setup_string += "Prune.fileName = 'skim_%s.root'\n" % (irfname)
    #else:
    #    setup_string += "Prune.fileName = 'skim_%s_%s.root'\n" % (
    #        irfname, livetime_bin)
    p8class = '%s_%s' % (class_version, classname)
    variant_cut = sel[variant]
    setup_string += "Prune.cuts = '(%s)&&(%s)'\n" % (variant_cut,
                                                     sel[classname])
    branches = get_branches(sel)
    branches += ['McEnergy', 'McXDir', 'McYDir', 'McZDir']
    branches += [config['Data']['var_xdir'], config['Data']['var_ydir'],
                 config['Data']['var_zdir'], config['Data']['var_energy']]
    setup_string += "Prune.branchNames = %s\n"%str(branches)

    merit_files = []    
    for k,v in generated.items():
        merit_files += sorted(glob.glob(os.path.join(config['meritdir'],
                                                     v['files'])))
    setup_string += 'Data.files = %s\n'%str(merit_files)
    setup_string += make_generated_events_str(generated)
    setup_string += """
Data.friends = {}
Data.tree_name = 'MeritTuple'\n
"""
    if livetime_bin is None:
        setup_string += "className='%s'\n" % (irfname)
    else:
        setup_string += "className='%s_%s'\n" % (irfname, livetime_bin)

    setup_string += "root_dir='.'\n"
    setup_string += write_config(config['Data'], 'Data')
    setup_string += write_config(config['Bins'], 'Bins')
    setup_string += write_config(config['FisheyeBins'], 'FisheyeBins')
    setup_string += write_config(config['EffectiveAreaBins'],
                                 'EffectiveAreaBins')
    setup_string += r"""
PSF.fit_pars={'ncore':[%g,1e-6,1.],
              'ntail':[%g,1e-6,1.0],
              'score':[%g,0.1,5.0],
              'stail':[%g,0.1,5],
              'gcore':[%g,1.001,20.],
              'gtail':[%g,1.001,20.]}
""" % tuple(config['PSF']['fit_pars'][variant])
    setup_string += ('PSF.scaling_pars=%s\n' %
                     str(config['PSF']['scaling_pars'][variant]))
    setup_string += ('Edisp.scaling_pars=%s\n' %
                     str(config['Edisp']['scaling_pars'][variant]))    
    if livetime_bin is None:
        setup_string += "parameterFile = 'par_%s.root'\n\n" % (irfname)
        setup_string += "selectionName=\'%s\'\n\n" % irfname
    else:
        setup_string += "parameterFile = 'par_%s_%s.root'\n\n" % (
            irfname, livetime_bin)
        setup_string += "selectionName=\'%s_%s\'\n\n" % (irfname, livetime_bin)


    setup_string += r"""
Bins.set_energy_bins()
Bins.set_angle_bins()
FisheyeBins.set_energy_bins()
FisheyeBins.set_angle_bins()
EffectiveAreaBins.set_energy_bins()
EffectiveAreaBins.set_angle_bins()
"""

    return setup_string


def main():

    usage = "%(prog)s [options] "
    description = """This generates directory structure and configuration
    files for IRF generation."""
    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('--meritdir',
                        default=None,
                        help='Set the root directory containing the input merit files.')

    parser.add_argument('--classdefs',
                        default=None, required=True,
                        help='Set the XML class definitions file.  This will be used to determine '
                        'the set of classes and event types for which IRFs will be generated.')

    parser.add_argument('--config',
                        default=None, required=True,
                        help='Set the configuration file.  This file defines scaling and fit '
                        'parameters as well as the energy range and number of thrown events in '
                        'each MC sample.')

    parser.add_argument('--livetime', default=False, action='store_true',
                        help='Peform setup for livetime correction.')

    args = parser.parse_args()

    # Default configuration
    config = {}
    selections = get_cuts_from_xml(args.classdefs)
    config.update(yaml.load(open(args.config, 'r')))

    if args.meritdir is not None:
        config['meritdir'] = args.meritdir

    # where to put the path structure
    target_dir = '.'
    class_version, classnames, typenames = get_class_types_from_xml(
        args.classdefs)

    config['MC'] = extract_generated_events(config['MC'], config['meritdir'])
    if 'Livetime' in config:
        config['Livetime']['MC'] = extract_generated_events(config['Livetime']['MC'],
                                                            config['meritdir'])
        
    
    for classname in classnames:
        for typename in typenames:
            irfname = '%s_%s_%s' % (class_version, classname,
                                    typename)
            print('Doing', irfname)
            dirname = os.path.join(target_dir, irfname)
            makedir(dirname)
            config['irf_name'] = str(irfname)
            config['irf_class'] = str(classname)
            config['irf_type'] = str(typename)
            with open(os.path.join(dirname,'config.yaml'),'w') as f:
                yaml.dump(config, f, Dumper=yaml.CDumper)

            with open(os.path.join(dirname, 'setup.py'), 'w') as f: 
                f.write(makesetup(class_version, classname,
                                  typename, irfname, config, selections))

        if not 'Livetime' in config:
            continue
                
        livetime_bins = sorted(config['Livetime']['MC'].keys())

        for (typename, ltbin) in itertools.product(typenames, livetime_bins):
            #irfname = '%s_%s_%s_%s' % (
            #    class_version, classname, irf_version, typename)
            irfname = '%s_%s_%s' % (class_version, classname, typename)
            print('Doing', irfname, ltbin)
            dirname = os.path.join(target_dir, irfname)
            f = open(os.path.join(dirname, 'setup_%s.py' % ltbin), 'w')
            f.write(makesetup(class_version, classname,
                              typename, irfname, config, selections, ltbin))
            f.close()


if __name__ == "__main__":
    main()
