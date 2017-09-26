"""
@file irfutils.py

@brief Utility functions for IRF generation.
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
