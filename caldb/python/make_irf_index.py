import os
import datetime
import numpy as np
from astropy.io import fits


def read_mappings(infile):
    """
    Read event_class or event_type mappings from a text file.
    Return matched sequences with the name to bit-position mapping.
    """
    event_class = [[], [], []]
    event_type = [[], [], []]
    for line in open(infile):
        # skip comment and whitespace-only lines
        if line[0] == "#" or not line.split():
            continue
        if line.startswith("event_class:"):
            mappings = event_class
            continue
        if line.startswith("event_type:"):
            mappings = event_type
            continue
        tokens = line.split()
        mappings[0].append(tokens[0])
        mappings[1].append(int(tokens[1]))
        try:
            mappings[2].append(int(tokens[2]))
        except ValueError:
            mappings[2].append(tokens[2])
    return event_class, event_type


def read_mappings_from_dict(index_mapping):
    """
    Read event_class and event_type mappings from a python dict.
    """
    evclass = [[], [], []]
    evtype = [[], [], []]
    for v in index_mapping['classes']:
        evclass[0].append(v[0])
        evclass[1].append(v[1])
        evclass[2].append(v[2])

    for v in index_mapping['types']:
        evtype[0].append(v[0])
        evtype[1].append(v[1])
        evtype[2].append(v[2])

    return evclass, evtype


def make_irf_index(index_mapping, outfile):

    if isinstance(index_mapping, dict):
        evclass, evtype = read_mappings_from_dict(index_mapping)
    else:
        evclass, evtype = read_mappings(index_mapping)

    if os.path.isfile(outfile):
        h_in = fits.open(outfile)

        for c in evclass[0]:
            m = h_in['BITMASK_MAPPING'].data.field('EVENT_CLASS') == c
            h_in['BITMASK_MAPPING'].data = h_in['BITMASK_MAPPING'].data[~m]

        evclass[0] = np.concatenate((h_in['BITMASK_MAPPING'].data.field('EVENT_CLASS'),
                                     evclass[0]))
        evclass[1] = np.concatenate((h_in['BITMASK_MAPPING'].data.field('BITPOSITION'),
                                     evclass[1]))
        evclass[2] = np.concatenate((h_in['BITMASK_MAPPING'].data.field('EVENT_TYPES'),
                                     evclass[2]))

    output = fits.HDUList()
    output.append(fits.PrimaryHDU())

    evclass_cols = [fits.Column('EVENT_CLASS', format='60A', unit=' ',
                                array=evclass[0]),
                    fits.Column('BITPOSITION', format='1I', unit=' ',
                                array=evclass[1]),
                    fits.Column('EVENT_TYPES', format='1J', unit='  ',
                                array=evclass[2])]
    evclass_hdu = fits.BinTableHDU.from_columns(
        evclass_cols, name='BITMASK_MAPPING')
    output.append(evclass_hdu)

    evtype_cols = [fits.Column('EVENT_TYPE', format='60A', unit=' ',
                               array=evtype[0]),
                   fits.Column('BITPOSITION', format='1I', unit=' ',
                               array=evtype[1]),
                   fits.Column('EVENT_TYPE_PARTITION', format='20A', unit=' ',
                               array=evtype[2])]
    evtype_hdu = fits.BinTableHDU.from_columns(
        evtype_cols, name='EVENT_TYPE_MAPPING')
    output.append(evtype_hdu)

    output.writeto(outfile, clobber=True)


if __name__ == '__main__':
    make_irf_index('irf_index.txt', 'irf_index_candidate.fits')
