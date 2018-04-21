
import mypy
import numpy
import pickle
from matplotlib import pyplot
import matplotlib.patheffects as PathEffects
from matplotlib.font_manager import FontProperties
from matplotlib.backends.backend_pdf import PdfPages

class combined_catalog(object):

	def __init__(self, header, n_galaxies):

		self.time_stamp = time.ctime()
		self.header = header

		self.n_galaxies = n_galaxies




catalog = pickle.load(open('../data/table_SFing_galaxies.pickle', 'rb'))






##########################################
###  Plot generic example UVJ diagram  ###
##########################################

if False:

    fig = pyplot.figure(figsize=(8.9, 8.))

    sp1 = fig.add_subplot(111)

    sp1.minorticks_on()

    sp1.set_xlabel('(V - J)$_{\mathrm{rest}}$', size=22)
    sp1.set_ylabel('(U - V)$_{\mathrm{rest}}$', size=22)

    sp1.axis([-0.2, 2.2, 0.2, 2.5])


    inds0 = (catalog.UVJ_class == 0) & (catalog.lmass > 8.25)
    inds1 = (catalog.UVJ_class == 1) & (catalog.lmass > 8.25)

    sp1.plot(catalog.VJ_color[inds0], catalog.UV_color[inds0], 'ko', ms=2, mec='r', zorder=1)
    sp1.plot(catalog.VJ_color[inds1], catalog.UV_color[inds1], 'ko', ms=2, mec='b', zorder=1)


    mypy.uvj_select_region(0.8, subplot=sp1, plot=True, kolor='#cccccc', lw=5,   zorder=2)
    mypy.uvj_select_region(0.8, subplot=sp1, plot=True, kolor='k', lw=2, zorder=3)


    font = FontProperties()
    font.set_family('sans-serif')


    t = sp1.text(0.03, 0.97, 'Quiescent', fontweight='bold', fontsize=22, color='w', 
    	         transform=sp1.transAxes, horizontalalignment='left', verticalalignment='top',
                 path_effects=[PathEffects.withStroke(linewidth=5, foreground='r')])

    t = sp1.text(0.88, 0.03, 'Star-\nForming', fontweight='bold', fontsize=22, color='w',  
    	         transform=sp1.transAxes, horizontalalignment='center', verticalalignment='bottom',
                 path_effects=[PathEffects.withStroke(linewidth=5, foreground='b')])








##################################################################
###  Plot UVJ distributions of each overdensity bin subsample  ###
##################################################################

if True:

    outname = '../figures/UVJ_overdens_bins.pdf'
    pdf = PdfPages(outname)



    fig = pyplot.figure(figsize=(8.9, 8.))

    sp1 = fig.add_subplot(111)

    sp1.minorticks_on()

    sp1.set_xlabel('(V - J)$_{\mathrm{rest}}$', size=22)
    sp1.set_ylabel('(U - V)$_{\mathrm{rest}}$', size=22)

    sp1.axis([-0.2, 2.2, 0.2, 2.5])

    mypy.uvj_select_region(0.8, subplot=sp1, plot=True, kolor='#cccccc', lw=5,   zorder=2)
    mypy.uvj_select_region(0.8, subplot=sp1, plot=True, kolor='k', lw=2, zorder=3)


    ###  text in upper-left corner
    text = 'All Q$_{spec}$>2 galaxies'
    text += '\nlog(M) > 8.5'
    t = sp1.text(0.03, 0.97, text, transform=sp1.transAxes,
                 horizontalalignment='left', verticalalignment='top')

    inds = (catalog.UVJ_class == 1) & (catalog.lmass > 8.5)
    points = sp1.plot(catalog.VJ_color[inds], catalog.UV_color[inds], 'ko', ms=2, mec='k', zorder=1)[0]

    pdf.savefig()


    colors = ['b', 'g', 'r']

    for i_densbin in range(3):

        text = catalog.overdens_labels[i_densbin]
        text += '\nlog(M) > 8.5'
        t.set_text(text)
        inds = (catalog.UVJ_class == 1) & (catalog.lmass > 8.5) & (catalog.digi_overdens == i_densbin+1)

        points.set_mec(colors[i_densbin])
        points.set_data((catalog.VJ_color[inds], catalog.UV_color[inds]))

        pdf.savefig()






    ###  iterating through all stellar mass bins
    for i_lmassbin in range(1, len(catalog.lmassbins)):

        text = 'All Q$_{spec}$>2 galaxies'
        text += '\nlog(M) = [%.2f-%.2f]' % (catalog.lmassbins[i_lmassbin-1], catalog.lmassbins[i_lmassbin])
        t.set_text(text)

        inds = (catalog.UVJ_class == 1) & (catalog.digi_lmass == i_lmassbin)
        points.set_mec('k')
        points.set_data((catalog.VJ_color[inds], catalog.UV_color[inds]))

        pdf.savefig()



        colors = ['b', 'g', 'r']

        for i_densbin in range(3):

            text = catalog.overdens_labels[i_densbin]
            text += '\nlog(M) = [%.2f-%.2f]' % (catalog.lmassbins[i_lmassbin-1], catalog.lmassbins[i_lmassbin])
            t.set_text(text)
            inds = (catalog.UVJ_class == 1) & (catalog.digi_lmass == i_lmassbin) & (catalog.digi_overdens == i_densbin+1)

            points.set_mec(colors[i_densbin])
            points.set_data((catalog.VJ_color[inds], catalog.UV_color[inds]))

            pdf.savefig()


    pdf.close()
    pyplot.close()
    print '\nwrote to:\n%s\n' % outname

























