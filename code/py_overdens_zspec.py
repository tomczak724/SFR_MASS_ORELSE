
import mypy
import numpy
import pickle
from scipy import optimize
from astropy.io import fits
from matplotlib import pyplot
import matplotlib.patheffects as PathEffects
from matplotlib.font_manager import FontProperties
from matplotlib.backends.backend_pdf import PdfPages




class combined_catalog(object):

    def __init__(self, header, n_galaxies):

        self.time_stamp = time.ctime()
        self.header = header

        self.n_galaxies = n_galaxies



catalog = pickle.load(open('../data/table_all_galaxies.pickle', 'rb'))













fig = pyplot.figure(figsize=(11.3, 6.45))

sp1 = fig.add_subplot(111)

sp1.minorticks_on()
sp1.set_xlabel('spectroscopic redshift')
sp1.set_ylabel('log( 1 + $\delta_{gal}$ )')

fig.subplots_adjust(left=0.11)

colors = ['b', 'g', 'r']



sp1.axis([0.66, 1.5, -0.4, 2.35])



sp1.plot(catalog.z, catalog.overdens, 'ko', mec='gray', ms=1)





'''
sp1.plot([0.7, 1.3], [0.5, 0.5], color='w', lw=4, ls='-', zorder=2)
sp1.plot([0.7, 1.3], [0.5, 0.5], color='k', lw=2, ls='-', zorder=2)

sp1.plot([0.7, 1.3], [1.0, 1.0], color='w', lw=4, ls='-', zorder=2)
sp1.plot([0.7, 1.3], [1.0, 1.0], color='k', lw=2, ls='-', zorder=2)

sp1.axvline(0.7, color='w', lw=4, ls='-', zorder=2)
sp1.axvline(0.7, color='k', lw=2, ls='-', zorder=2)

sp1.axvline(1.3, color='w', lw=4, ls='-', zorder=2)
sp1.axvline(1.3, color='k', lw=2, ls='-', zorder=2)
'''




sp1.axhline(0.5, color='w', lw=4, ls='-', zorder=2)
sp1.axhline(0.5, color='k', lw=2, ls='-', zorder=2)

sp1.axhline(1.0, color='w', lw=4, ls='-', zorder=2)
sp1.axhline(1.0, color='k', lw=2, ls='-', zorder=2)






labels = ['field', 'groups', 'clusters']

for i_overdens in range(len(catalog.overdensbars)):

    inds = (catalog.digi_overdens == i_overdens+1)


    if i_overdens in [0, 2]:
        x_text, y_text = sp1.transLimits.transform((1.4, numpy.median(catalog.overdens[inds])))
    else:
        x_text, y_text = sp1.transLimits.transform((1.4, 0.75))

    sp1.text(x_text, y_text, '"%s"\nN = %i' % (labels[i_overdens], numpy.count_nonzero(inds)), 
             transform=sp1.transAxes, horizontalalignment='center', verticalalignment='center', 
             color=colors[i_overdens], fontweight='normal', fontsize=20, 
             path_effects=[PathEffects.withStroke(linewidth=1, foreground='k')])









