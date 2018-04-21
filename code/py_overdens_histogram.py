
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













fig = pyplot.figure(figsize=(7.2, 6.5))

sp1 = fig.add_subplot(111)

sp1.minorticks_on()
sp1.set_xlabel('log( 1 + $\delta_{gal}$ )')
sp1.set_ylabel('Number')



sp1.set_yscale('log')
sp1.axis([-0.8, 2.3, 0.2, 2000])
sp1.set_yticklabels(['', '', '1 ', '10 ', '100 ', '1000 '])

fig.subplots_adjust(left=0.14, bottom=0.11, right=0.99, top=0.99)




nbins = 58
h0 = sp1.hist(catalog.overdens, bins=nbins, range=(-0.7, 2.2), 
              histtype='step', color='k', lw=2, log=0, zorder=2)


#sp1.axvline(0.5, color='k', lw=2)
#sp1.axvline(1.0, color='k', lw=2)





colors = ['b', 'g', 'r']
colors = ['#3f4ea1', '#6db388', '#d92120']

for i_overdens in range(len(catalog.overdensbars)):

    inds = (catalog.digi_overdens == i_overdens+1)

    h1 = sp1.hist(catalog.overdens[inds], bins=nbins, range=(-0.7, 2.2), lw=2,
                  histtype='stepfilled', color=colors[i_overdens], log=0, zorder=2)







t1 = sp1.text(-0.37, 170, 'low\ndensity', fontsize=20,
              horizontalalignment='center', verticalalignment='bottom') 
t2 = sp1.text(0.85, 250, 'intermediate\ndensity', fontsize=20,
              horizontalalignment='center', verticalalignment='bottom') 
t3 = sp1.text(1.75, 50, 'high\ndensity', fontsize=20,
              horizontalalignment='center', verticalalignment='bottom')


















