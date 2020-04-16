"""
@file irfutils.py

@brief Utility functions for IRF generation.
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import yaml
import re
import collections
import glob
import array
from xml.dom import minidom
import xml.etree.cElementTree as ElementTree
import argparse
import itertools
import numpy as np
import ROOT
from astropy.io import fits


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
            cut_dict[k] = cut_dict[k].replace(k0, '(%s)' % v0)


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


class read_aeff(object):
    def __init__(self, fname):
        hdulist = fits.open(fname)
        self.hdu1 = hdulist[1].data
        self.e_lo = self.hdu1.field('ENERG_LO')[0]
        self.e_hi = self.hdu1.field('ENERG_HI')[0]
        self.cth_lo = self.hdu1.field('CTHETA_LO')[0]
        self.cth_hi = self.hdu1.field('CTHETA_HI')[0]
        self.aeff = self.hdu1.field('EFFAREA')[0]
        self.n_e = len(self.e_lo)
        self.n_cth = len(self.cth_lo)
        print('AEFF TABLE', self.n_e, '*', self.n_cth)

    def get_e_bin_bounds(self, i):
        if i < 0 or i >= self.n_e:
            return -1
        return (self.e_lo[i], self.e_hi[i])

    def get_loge_bin_bounds(self, i):
        if i < 0 or i >= self.n_e:
            return -1
        return (np.log10(self.e_lo[i]), np.log10(self.e_hi[i]))

    def get_e_bin_center(self, i):
        if i < 0 or i >= self.n_e:
            return -1
        return (self.e_lo[i] + self.e_hi[i]) / 2

    def get_a_bin_bounds(self, i):
        if i < 0 or i >= self.n_cth:
            return -1
        return [degrees(acos(self.cth_hi[i])),
                degrees(acos(self.cth_lo[i]))]

    def get_a_bin_center(self, i):
        if i < 0 or i >= self.n_cth:
            return -1
        return degrees(acos((self.cth_lo[i] + self.cth_hi[i]) / 2))

    def get_aeff(self, ebin, abin):
        return self.aeff[abin, ebin]

    def get_averaged_aeff(self, ebin):
        a = 0
        for i in range(self.n_cth):
            a += self.get_aeff(ebin, i)
        a /= self.n_cth
        return a

    def get_closest_e_bin(self, e):
        for i in range(self.n_e):
            if self.e_lo[i] <= e and self.e_hi[i] > e:
                return i
        return -1


class lockedfit(object):
    def __init__(self, avltime):
        self.x0 = avltime

    def __call__(self, x, par):
        m = par[0]
        x1 = x[0] - self.x0
        return m * x1 + 1

    def intercept(self, m):
        return -m * self.x0 + 1


class p01fit(object):
    def __call__(self, var, par):
        (a0, b0, a1, logEb1, a2, logEb2) = (
            par[0], par[1], par[2], par[3], par[4], par[5])
        logE = var[0]

        b1 = (a0 - a1) * logEb1 + b0
        b2 = (a1 - a2) * logEb2 + b1
        if logE < logEb1:
            return a0 * logE + b0
        if logE < logEb2:
            return a1 * logE + b1
        return a2 * logE + b2


class corr_fit(object):
    def __init__(self, base_aeff_fname, av_ltime, canvas_base_name="c0",
                 aeff_fnames_and_ltimes=[], min_ebin=10, max_ebin=70):
        self.base_aeff = read_aeff(base_aeff_fname)
        self.aeff_set = []
        self.ltimes = []
        self.cname = canvas_base_name
        if len(aeff_fnames_and_ltimes) > 0:
            add_aeff_files_array(aeff_fnames_and_ltimes)
        # default energy bins: (10) ~ 100 MeV to (max) ~ 500 GeV
        self.ebins = []
        self.av_ltime = av_ltime
        for i in range(max_ebin - min_ebin):
            self.ebins.append(min_ebin + i)
        self.n_ebins = len(self.ebins)
        # setup sim evt numbers: n events from lemin to lemax
        self.use_nsims = True
        self.nsims = []
        self.sim_lemin = []
        self.sim_lemax = []
        # any other setup to do
        ROOT.gStyle.SetOptFit(111) #was 0111 but in python 3 that is an octal number 
        ROOT.gROOT.SetStyle('Plain')

    def add_aeff_file(self, aeff_fname_and_ltime, sim_ev_data=[]):
        fname = aeff_fname_and_ltime[0]
        ltime = aeff_fname_and_ltime[1]
        if len(sim_ev_data) == 0:
            self.use_nsims = False

        if self.use_nsims:
            self.nsims.append(sim_ev_data[0])
            self.sim_lemin.append(sim_ev_data[1])
            self.sim_lemax.append(sim_ev_data[2])

        self.aeff_set.append(read_aeff(fname))
        self.ltimes.append(ltime)

    def add_aeff_files_array(self, ae_array, sim_ev_array=[]):
        for l in range(len(ae_array)):
            if len(sim_ev_array) > 0:
                add_aeff_file(ae_array[l], sim_ev_array[l])
            else:
                add_aeff_file(ae_array[l], [])

    def n_ltbins(self):
        return len(self.ltimes)

    def fit(self,
            old_p0_pars=[0., 0., 0., 0., 0., 0.],
            old_p0_plimits=[],
            fits_output=""):
        #
        # for each energy bin fits angle-averaged aeff
        # versus livetime bin center with a pol1
        # Linear fit parameters are plotted versus energy
        # and fitted with a piecewise linear (3 pieces)
        #

        self.c1 = ROOT.TCanvas(self.cname + '_c1', self.cname, 900, 800)
        self.c1.Divide(7, 10)

        self.c3 = []
        for i in range(4):
            c = ROOT.TCanvas(self.cname + '_c%i' % (i + 3),
                             self.cname + '_c%i' % (i), 900, 800)
            c.Divide(4, 4)
            self.c3.append(c)

        # set up global histograms
        self.g_p0 = ROOT.TGraphErrors(7)
        self.g_p1 = ROOT.TGraphErrors(7)
        self.graphs = []
        if len(fits_output) > 0:
            fits_e = []
            fits_p0 = []
            fits_p1 = []

        # do a aeff_vs_ltime plot and fit, for each energy bin
        for i in range(self.n_ebins):
            self.c1.cd(i + 1)
            # if no sim data infos: simple tgraph
            if not self.use_nsims:
                g = ROOT.TGraph(self.n_ltbins())
            else:
                g = ROOT.TGraphErrors(self.n_ltbins())
            bin = self.ebins[i]
            base_aeff = self.base_aeff.get_averaged_aeff(self.ebins[i])
            for li in range(self.n_ltbins()):
                this = self.aeff_set[li].get_averaged_aeff(self.ebins[i])
                g.SetPoint(li, self.ltimes[li], this / base_aeff)
                # if sim data info esist, setup errors accordingly
                if self.use_nsims:
                    # loop across sim intervals, calculate sim evts in interval
                    nevts = 0.
                    for si in range(len(self.nsims)):
                        (bin_lemin, bin_lemax) = self.aeff_set[li].get_loge_bin_bounds(
                            i)
                        nevts += 1. * \
                            self.nsims[si] / (self.sim_lemax[si] -
                                              self.sim_lemin[si]) * (bin_lemax - bin_lemin)
                    ndetected = nevts * this / 6.  # ndet=nsim*efficiency
                    err = np.sqrt(ndetected) / ndetected * this / base_aeff
                    g.SetPointError(li, 0, err)
                    # if i==2 or i==40:
                    # print 'point',i,':',self.ltimes[li],'ndet,nsim',ndetected,nevts,'aeff,ratio',this,this/base_aeff,'sqrt(n)/n,err',sqrt(ndetected)/ndetected,err,'E',self.base_aeff.get_e_bin_center(self.ebins[i])
            bounds = self.base_aeff.get_e_bin_bounds(self.ebins[i])
            center = self.base_aeff.get_e_bin_center(self.ebins[i])
            name = "Energy: [%.1f,%.1f](%.1f)" % (bounds[0], bounds[1], center)
            g.SetTitle(name)
            g.Draw("AP")
            g.SetMarkerStyle(7)
            # linear fit, locked to have 1 at average livetime
            lffunc = lockedfit(self.av_ltime)
            lf1 = ROOT.TF1("lffunc", lffunc, -1, 1, 1)
            lf1.SetParameter(0, 1)
            g.Fit(lf1)
            func = g.GetFunction('lffunc')
            func.SetLineColor(2)
            func.SetLineWidth(1)
            # store fit parameters in global graphs
            fitted_m = func.GetParameter(0)
            fitted_q = lffunc.intercept(fitted_m)
            print('>>>>>', fitted_m, fitted_q)
            self.g_p0.SetPoint(i, np.log10(center), fitted_m)
            self.g_p1.SetPoint(i, np.log10(center), fitted_q)
            w = 0  # (np.log10(bounds[1])-np.log10(bounds[0]))/np.sqrt(12)
            error_m = func.GetParError(0)
            # error_q=abs(error_m/fitted_m*fitted_q)
            error_q = abs(error_m * fitted_q / fitted_m * self.av_ltime)
            self.g_p0.SetPointError(i, w, error_m)
            self.g_p1.SetPointError(i, w, error_q)
            self.graphs.append(g)
            if len(fits_output) > 0:
                fits_e.append(center)
                fits_p0.append(func.GetParameter(1))
                fits_p1.append(func.GetParameter(0))

            # Draw on c3 canvas
            j = i // 16

            if j >= len(self.c3):
                continue

            self.c3[j].cd(i % 16 + 1)
            g.Draw("AP")
            func.Draw("SAME")

        self.c1.Update()
        self.c1.SaveAs('livetime_fits_' + self.cname + '.png')

        for i, c in enumerate(self.c3):
            c.Update()
            c.SaveAs('livetime_fits%i_' % (i) + self.cname + '.png')

        # save results in Fits file
        if len(fits_output) > 0:
            col1 = fits.Column(name='ENERGY', format='E',
                               array=numpy.array(fits_e))
            col2 = fits.Column(name='P0', format='E',
                               array=numpy.array(fits_p0))
            col3 = fits.Column(name='P1', format='E',
                               array=numpy.array(fits_p1))
            cols = fits.ColDefs([col1, col2, col3])
            hdu2 = fits.new_table(cols)
            hdu1 = fits.PrimaryHDU()
            hdulist = fits.HDUList([hdu1, hdu2])
            hdulist.writeto(fits_output, clobber=True)

        # canvas for global fit
        self.c2 = ROOT.TCanvas(self.cname + '_c2', self.cname, 900, 450)
        self.c2.Divide(2, 1)
        self.c2.cd(1)
        self.g_p0.SetTitle("P0")
        self.g_p0.Draw("AP")
        self.g_p0.SetMarkerStyle(7)
        self.g_p0.GetXaxis().SetTitle('log(E)')

        #### P0 ####
        f_p0 = p01fit()
        self.tf1_p0 = ROOT.TF1('tf1_p0', f_p0, 2, 5.5, 6)
        self.tf1_p0.SetParNames('a0', 'b0', 'a1', 'logEb1', 'a2', 'logEb2')
        # start from old vals
        self.tf1_p0.SetParameter(0, old_p0_pars[0])  # a0
        self.tf1_p0.SetParameter(1, old_p0_pars[1])  # b0
        self.tf1_p0.SetParameter(2, old_p0_pars[2])  # a1
        self.tf1_p0.SetParameter(3, old_p0_pars[3])  # Eb1
        self.tf1_p0.SetParameter(4, old_p0_pars[4])  # a2
        self.tf1_p0.SetParameter(5, old_p0_pars[5])  # Eb2
        self.tf1_p0.SetParLimits(3, 1.5, 2.5)
        self.tf1_p0.SetParLimits(5, 2.5, 4.0)
        if len(old_p0_plimits) > 0:
            for plmi in range(len(old_p0_plimits)):
                self.tf1_p0.SetParLimits(
                    plmi, old_p0_plimits[plmi][0], old_p0_plimits[plmi][1])
        self.g_p0.Fit(self.tf1_p0)
        self.tf1_p0.SetLineWidth(1)
        self.tf1_p0.SetLineColor(ROOT.kRed)
        self.c2.Update()

        self.c2.cd(2)
        self.g_p1.SetTitle("P1")
        self.g_p1.Draw("AP")
        self.g_p1.SetMarkerStyle(7)
        self.g_p1.GetXaxis().SetTitle('log(E)')
        #### P1 ####
        f_p1 = p01fit()
        self.tf1_p1 = ROOT.TF1('tf1_p1', f_p1, 1.5, 5.5, 6)
        self.tf1_p1.SetParNames('a0', 'b0', 'a1', 'logEb1', 'a2', 'logEb2')
        # old vals
        new_a0 = self.tf1_p0.GetParameter(0)
        new_b0 = self.tf1_p0.GetParameter(1)
        new_a1 = self.tf1_p0.GetParameter(2)
        new_eb1 = self.tf1_p0.GetParameter(3)
        new_a2 = self.tf1_p0.GetParameter(4)
        new_eb2 = self.tf1_p0.GetParameter(5)
        self.final_results = [[new_a0, new_b0, new_a1,
                               new_eb1, new_a2, new_eb2]]
        new_a0 = -new_a0 * self.av_ltime
        new_b0 = 1 - new_b0 * self.av_ltime
        new_a1 = -new_a1 * self.av_ltime
        new_a2 = -new_a2 * self.av_ltime
        self.final_results.append([new_a0, new_b0, new_a1,
                                   new_eb1, new_a2, new_eb2])

        self.tf1_p1.SetParameter(0, new_a0)  # a0
        self.tf1_p1.SetParameter(1, new_b0)  # b0
        self.tf1_p1.SetParameter(2, new_a1)  # a1
        self.tf1_p1.SetParameter(3, new_eb1)  # Eb1
        self.tf1_p1.SetParameter(4, new_a2)  # a2
        self.tf1_p1.SetParameter(5, new_eb2)  # Eb2
        self.tf1_p1.Draw("SAME")
        self.tf1_p1.SetLineWidth(1)
        self.tf1_p1.SetLineColor(ROOT.kBlue)

        self.g_p1.GetHistogram().SetMaximum(3.0)
        self.g_p1.GetHistogram().SetMinimum(-3.0)

        self.c2.Update()
        self.c2.SaveAs('livetime_params_' + self.cname + '.png')

        print('P1  A0', new_a0)
        print('P1  B0', new_b0)
        print('P1  A1', new_a1)
        print('P1 EB1', new_eb1)
        print('P1  A2', new_a2)
        print('P1 EB2', new_eb2)
        print('done')

    def plot_angle_dep(self, angle_bins=[3, 7, 11, 15, 19, 23, 27, 31]):
        #
        # for selected angle bins do the linear fit of
        # aeff versus livetime, plots the pol1 parameters
        # versus energy in a plot to evaluate angle
        # dependence of correction
        #
        self.fitcanvases = []
        self.angle_graphs = []
        self.angle_globals_p0 = []
        self.angle_globals_p1 = []
        for ci in range(len(angle_bins)):
            a = self.base_aeff.get_a_bin_center(angle_bins[ci])
            cname = "%s_c3_%d" % (self.cname, ci)
            ctit = "angle: %.1f" % (a)
            print('doing', ctit)
            canv = ROOT.TCanvas(cname, ctit, 900, 800)
            canv.Divide(8, 7)
            # set up global histograms
            g_a_p0 = ROOT.TGraphErrors(7)
            g_a_p1 = ROOT.TGraphErrors(7)
            # do a aeff_vs_ltime plot and fit, for each energy bin
            for i in range(self.n_ebins):
                canv.cd(i + 1)
                g = ROOT.TGraph(self.n_ltbins())
                base_aeff = self.base_aeff.get_aeff(
                    self.ebins[i], angle_bins[ci])
                for li in range(self.n_ltbins()):
                    this = self.aeff_set[li].get_aeff(
                        self.ebins[i], angle_bins[ci])
                    g.SetPoint(li, self.ltimes[li], this / base_aeff)
                bounds = self.base_aeff.get_e_bin_bounds(self.ebins[i])
                center = self.base_aeff.get_e_bin_center(self.ebins[i])
                name = "Energy: [%.1f,%.1f](%.1f)" % (
                    bounds[0], bounds[1], center)
                g.SetTitle(name)
                g.Draw("AP")
                g.SetMarkerStyle(7)
                # linear fit
                g.Fit('pol1')
                func = g.GetFunction('pol1')
                func.SetLineColor(2)
                func.SetLineWidth(1)
                # store fit parameters in global graphs
                g_a_p1.SetPoint(i, np.log10(center), func.GetParameter(0))
                g_a_p0.SetPoint(i, np.log10(center), func.GetParameter(1))
                w = 0  # (np.log10(bounds[1])-np.log10(bounds[0]))/np.sqrt(12)
                g_a_p1.SetPointError(i, w, func.GetParError(0))
                g_a_p0.SetPointError(i, w, func.GetParError(1))
                # persistency
                self.angle_graphs.append(g)
            self.fitcanvases.append(canv)
            # store globals
            g_a_p0.SetTitle('P0')
            g_a_p1.SetTitle('P1')
            self.angle_globals_p0.append(g_a_p0)
            self.angle_globals_p1.append(g_a_p1)
        # canvas for global plots
        self.c4 = ROOT.TCanvas(self.cname + '_c4', self.cname, 900, 450)
        self.c4.Divide(2, 1)
        self.c4.cd(1)  # P0#####
        self.leg0_a = ROOT.TLegend(.65, .65, .95, .95)
        for gi in range(len(self.angle_globals_p0)):
            if gi == 0:
                plotop = "AP"
            else:
                plotop = "P"
            self.angle_globals_p0[gi].Draw(plotop)
            self.angle_globals_p0[gi].SetLineColor(gi + 1)
            tit = "angle: %.1f" % (
                self.base_aeff.get_a_bin_center(angle_bins[gi]))
            self.leg0_a.AddEntry(self.angle_globals_p0[gi], tit, "lp")
        self.leg0_a.Draw()

        self.c4.cd(2)  # P1#####
        self.leg1_a = ROOT.TLegend(.65, .15, .95, .45)
        for gi in range(len(self.angle_globals_p1)):
            if gi == 0:
                plotop = "AP"
            else:
                plotop = "P"
            self.angle_globals_p1[gi].Draw(plotop)
            self.angle_globals_p1[gi].SetLineColor(gi + 1)
            tit = "angle: %.1f" % (
                self.base_aeff.get_a_bin_center(angle_bins[gi]))
            self.leg1_a.AddEntry(self.angle_globals_p1[gi], tit, "lp")
        self.leg1_a.Draw()

    def writein(self, fitsname):
        f = fits.open(fitsname)
        newcol = fits.Column(name='EFFICIENCY_PARS', format='6E',
                             array=self.final_results)
        newdata = fits.FITS_rec.from_columns([newcol])
        table_hdu = fits.BinTableHDU(newdata, name='EFFICIENCY_PARAMS')
        headkeys = []
        for key in ['EXTNAME', 'TELESCOP', 'INSTRUME', 'FILTER',
                    'HDUCLASS', 'HDUCLAS1', 'HDUCLAS2', 'HDUVERS',
                    'CCLS0001', 'CDTP0001', 'CCNM0001', 'CBD10001',
                    'CBD20001', 'CBD30001', 'CBD40001', 'CBD50001',
                    'CBD60001', 'CBD70001', 'CBD80001', 'CBD90001',
                    'CVSD0001', 'CVST0001', 'CDES0001']:
            headkeys.append([key, f[3].header[key]])
        f[3] = table_hdu
        for pars in headkeys:
            f[3].header[pars[0]] = pars[1]

        f.writeto(fitsname, clobber=True, checksum=True)
        f.close()
