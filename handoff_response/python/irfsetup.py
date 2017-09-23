"""
@file irfsetup.py

@brief Generate python setup files and directory structure for
makeirf/makefits.
"""

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


def get_branches(aliases):
    """Get unique branch names from an alias dictionary."""
    
    ignore = ['pow', 'log10', 'sqrt', 'max']
    branches = []
    for k, v in aliases.items():

        tokens = re.sub('[\(\)\+\*\/\,\=\<\>\&\!\-\|]', ' ', v).split()

        for t in tokens:

            if bool(re.search(r'^\d', t)) or len(t) <= 3:
                continue

            if bool(re.search(r'[a-zA-Z]', t)) and t not in ignore:
                branches += [t]

    return list(set(branches))


def update_dict(d, u):
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            r = update_dict(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d


def write_config(config, name):

    out_str = ''
    for k, v in sorted(config.items()):
        if isinstance(v, str):
            out_str += r"{}.{} = '{}'".format(name, k, v) + '\n'
        else:
            out_str += r"{}.{} = {}".format(name, k, v) + '\n'
    return out_str


def extract_generated_events(config, meritdir):
    """Extract the number of generated events in each MC sample.  If
    'ngen' is not defined this method will try to extract this
    information from the 'jobinfo' tree.

    Parameters
    ----------
    config : dict
        Configuration dictionary containing one dictionary for each MC
        sample.
    """

    cfg_out = {}

    for k, v in config.items():

        if isinstance(v, dict):

            cfg_out[k] = v.copy()

            if 'ngen' in v:
                continue
            elif 'files' in v:
                chain = ROOT.TChain('jobinfo')
                for f in glob.glob(os.path.join(meritdir, v['files'])):
                    chain.Add(f)
                cfg_out[k]['ngen'] = getGeneratedEvents(chain)
            else:
                raise ValueError(r"Either 'ngen' or 'files' must be defined.")
        elif isinstance(v, list):

            ngen = eval(v[0])
            logemin = eval(str(v[1]))
            logemax = eval(str(v[2]))
            cfg_out[k] = {'ngen': ngen, 'logemin': logemin, 'logemax': logemax}

        else:
            raise ValueError('Invalid type for MC dictionary element.')

    return cfg_out


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


def getGeneratedEvents(chain):
    NGen_sum = 0
    vref = {}

    vref['trigger'] = array.array('i', [0])
    vref['generated'] = array.array('i', [0])
    vref['version'] = array.array('f', [0.0])
    vref['revision'] = array.array('f', [0.0])
    vref['patch'] = array.array('f', [0.0])
    chain.SetBranchAddress('trigger', vref['trigger'])
    chain.SetBranchAddress('generated', vref['generated'])
    chain.SetBranchAddress('version', vref['version'])

    if chain.GetListOfBranches().Contains('revision'):
        chain.SetBranchAddress('revision', vref['revision'])

    if chain.GetListOfBranches().Contains('patch'):
        chain.SetBranchAddress('patch', vref['patch'])

    for i in range(chain.GetEntries()):
        chain.GetEntry(i)
        ver = int(vref['version'][0])
        rev = int(vref['revision'][0])
        patch = int(vref['patch'][0])
        NGen = vref['generated'][0]
        NGen_sum += NGen

    return NGen_sum


def strip(input_str):
    return str(input_str.replace('\n', '').replace(' ', ''))


def replace_aliases(cut_dict, aliases):

    for k, v in cut_dict.items():
        for k0, v0 in aliases.items():
            cut_dict[k] = cut_dict[k].replace(k0, v0)


def makedir(dirname):

    if os.path.exists(dirname):
        return

    try:
        os.mkdir(dirname)
    except OSError as err:
        if update:
            pass  # do not care if direcotry exists
        else:
            raise err


def get_class_types_from_xml(xmlfile):
    xmldoc = minidom.parse(xmlfile)
    class_version = xmldoc.getElementsByTagName(
        'EventClass')[0].attributes['version'].value
    itemlist = xmldoc.getElementsByTagName('EventMap')
    classnames = [l.attributes['name'].value for l in itemlist[0].getElementsByTagName(
        "EventCategory") if l.attributes['name'].value != "LLE"]
    typenames = [l.attributes['name'].value for l in itemlist[1].getElementsByTagName(
        "EventCategory")]
    return class_version, classnames, typenames


def get_cuts_from_xml(xmlfile):
    """Extract event selection strings from the XML file."""

    root = ElementTree.ElementTree(file=xmlfile).getroot()
    event_maps = root.findall('EventMap')
    alias_maps = root.findall('AliasDict')[0]

    event_classes = {}
    event_types = {}
    event_aliases = {}

    for m in event_maps:
        if m.attrib['altName'] == 'EVENT_CLASS':
            for c in m.findall('EventCategory'):
                event_classes[c.attrib['name']] = strip(
                    c.find('ShortCut').text)
        elif m.attrib['altName'] == 'EVENT_TYPE':
            for c in m.findall('EventCategory'):
                event_types[c.attrib['name']] = strip(c.find('ShortCut').text)

    for m in alias_maps.findall('Alias'):
        event_aliases[m.attrib['name']] = strip(m.text)

    replace_aliases(event_aliases, event_aliases.copy())
    replace_aliases(event_aliases, event_aliases.copy())
    replace_aliases(event_classes, event_aliases)
    replace_aliases(event_types, event_aliases)

    event_selections = {}
    event_selections.update(event_classes)
    event_selections.update(event_types)
    event_selections.update(event_aliases)

    return event_selections


def makesetup(class_version, classname, variant, irfname, config, livetime_bin=None):

    sel = config['selections']
    generated = config['MC']

    if not livetime_bin:
        globfile = '*.root'
        meritdir = config['meritdir']
    else:
        meritdir = os.path.join(config['meritdir'], 'livetimebins')
        globfile = '*%s*.root' % livetime_bin

    setup_string = """
from IRFdefault import *
import glob
from math import *\n
"""

    setup_string += "Prune.fileName = ''"    
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
    setup_string += """
meritFiles = '%s'
Data.files = sorted(glob.glob(meritFiles))\n
""" % (os.path.join(meritdir, globfile))

    if livetime_bin is None:
        setup_string += make_generated_events_str(generated)
    else:
        n_generated = float(config['Livetime']['ngen'][livetime_bin])
        setup_string += "Data.generated = [%g]\n" % n_generated
        setup_string += "Data.logemin = [1.25]\n"
        setup_string += "Data.logemax = [5.75]\n"

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

    if livetime_bin is None:
        setup_string += "parameterFile = 'par_%s.root'\n\n" % (irfname)
        setup_string += "selectionName=\'%s\'\n\n" % irfname
    else:
        setup_string += "parameterFile = 'par_%s_%s.root'\n\n" % (
            irfname, livetime_bin)
        setup_string += "selectionName=\'%s_%s\'\n\n" % (irfname, livetime_bin)
    setup_string += ('Edisp.scaling_pars=%s\n' %
                     str(config['Edisp']['scaling_pars'][variant]))

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

    default_bins = {'logemin': 0.75, 'logemax': 6.5, 'logedelta': 0.25,
                    'cthmin': 0.2, 'cthmax': 1.0, 'cthdelta': 0.1,
                    'edisp_energy_overlap': 0, 'edisp_angle_overlap': 0,
                    'psf_energy_overlap': 0, 'psf_angle_overlap': 0}

    # Default configuration
    config = {}
    config['selections'] = get_cuts_from_xml(args.classdefs)
    config.update(yaml.load(open(args.config, 'r')))

    if args.meritdir is not None:
        config['meritdir'] = args.meritdir

    # where to put the path structure
    target_dir = '.'
    class_version, classnames, typenames = get_class_types_from_xml(
        args.classdefs)

    irf_version = config['irf_version']
    config['MC'] = extract_generated_events(config['MC'], config['meritdir'])

    for classname in classnames:
        for typename in typenames:
            irfname = '%s_%s_%s_%s' % (class_version, classname,
                                       irf_version, typename)
            print 'Doing', irfname
            dirname = target_dir + '/' + irfname + '/'
            makedir(dirname)
            f = open(dirname + 'setup.py', 'w')
            f.write(makesetup(class_version, classname, typename, irfname, config))
            f.close()

    if args.livetime:
        livetime_bins = config['Livetime']['ngen'].keys()

        for classname in classnames:
            for (typename, ltbin) in itertools.product(typenames, livetime_bins):
                irfname = '%s_%s_%s_%s' % (
                    class_version, classname, irf_version, typename)
                print 'Doing', irfname, ltbin
                dirname = target_dir + '/' + irfname + '/'
                makedir(dirname)
                f = open(dirname + 'setup_%s.py' % ltbin, 'w')
                f.write(makesetup(class_version, classname,
                                  typename, irfname, config, ltbin))
                f.close()


if __name__ == "__main__":
    main()
