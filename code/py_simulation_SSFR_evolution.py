
import time
import fsps
import mypy
import glob
import numpy
import progressbar
from os import path
from astropy import cosmology
from matplotlib import pyplot, patheffects
from matplotlib.backends.backend_pdf import PdfPages

cosmo = cosmology.FlatLambdaCDM(H0=70., Om0=0.3)

###  establishing redshift-age relation for this cosmology
print '\nGenerating cosmology...'
t0 = time.time()
universe_redshifts = numpy.linspace(30, 0, 2000)
universe_age = cosmo.age(universe_redshifts)
tf = time.time()
print 'done!  %.1f seconds\n' % (tf-t0)



'''
###  setting up FSPS model
print '\nGenerating FSPS stellar population model...'
t0 = time.time()

sps = fsps.StellarPopulation(zcontinuous=1)  # zcontinuous = 1 = interpolate SSPs to value given by "logzsol"
sps.params['logzsol'] = 0                 # Adopting solar metallicity
sps.params['sfh'] = 3                     # 1 = Tau model ; 3 = User input model
sps.params['imf_type'] = 1                # Chabrier (2003)
sps.params['add_igm_absorption'] = True   # Add IGM absorption
sps.params['dust_type'] = 2               # Calzetti+2000 attenuation curve

filternames = ['U', 'V', '2MASS_J']

tf = time.time()
print 'done!  %.1f seconds\n' % (tf-t0)
'''




def log_SFR_tomczak2016(z, lmass, SF_galaxies=True):
    '''
    Star-forming galaxies:
      s0      = 0.448 + 1.220*z - 0.174*z**2
      log(M0) = 9.458 + 0.865*z - 0.132*z**2
      gamma   = 1.091

    All galaxies:
      s0      = 0.195 + 1.157*z - 0.143*z**2
      log(M0) = 9.244 + 0.753*z - 0.090*z**2
      gamma   = 1.118

    '''

    if SF_galaxies:
        s0_params = [0.448, 1.220, -0.174]
        m0_params = [9.458, 0.865, -0.132]
        gamma = 1.091
    else:
        s0_params = [0.195, 1.157, -0.143]
        m0_params = [9.244, 0.753, -0.090]
        gamma = 1.118

    s0 = s0_params[0] + s0_params[1] * z + s0_params[2] * z**2
    log_m0 = m0_params[0] + m0_params[1] * z + m0_params[2] * z**2
    return s0 - numpy.log10(1 + (10**lmass / 10**log_m0)**-gamma)


def dschechter(lmassax, lmstar, a1, a2, phistar1, phistar2):
        factor1 = numpy.log(10) * numpy.exp(-10**(lmassax - lmstar)) * 10**(lmassax - lmstar)
        factor2 = phistar1 * 10**(a1*(lmassax - lmstar)) + phistar2 * 10**(a2*(lmassax - lmstar))
        return factor1 * factor2


def fquiescent_leja2015(z, lmass):
        '''
        Returns the Quiescent fraction at the given redshift
        and stellar mass making use of the parameterizations
        from Leja+2015. In brief, this function calculates
        the QUIESCENT and STAR-FORMING stellar mass functions
        at z, lmass and returns the ratio: QU / (QU + SF)
        '''
        #logphi1_tot  = -2.64 + 0.07 * z - 0.28 * z**2
        #logphi2_tot  = -3.11 - 0.18 * z - 0.03 * z**2
        #logMstar_tot = 10.72 - 0.13 * z + 0.11 * z**2
        #alpha1_tot = -0.39
        #alpha2_tot = -1.53

        logphi1_sf  = -2.88 + 0.11 * z - 0.31 * z**2
        logphi2_sf  = -3.48 + 0.07 * z - 0.11 * z**2
        logMstar_sf = 10.67 - 0.02 * z + 0.10 * z**2
        alpha1_sf = -0.97
        alpha2_sf = -1.58

        logphi1_qu  = -2.51 - 0.33 * z - 0.07 * z**2
        logphi2_qu  = -3.54 - 2.31 * z + 0.27 * z**2
        logMstar_qu = 10.70
        alpha1_qu = -0.10
        alpha2_qu = -1.69

        smf_sf = dschechter(lmass, logMstar_sf, alpha1_sf, alpha2_sf, 10**logphi1_sf, 10**logphi2_sf)
        smf_qu = dschechter(lmass, logMstar_qu, alpha1_qu, alpha2_qu, 10**logphi1_qu, 10**logphi2_qu)

        return smf_qu / (smf_qu + smf_sf)








###  setting global variables
DT = 0.01
GRID_TIME = numpy.arange(DT, 12., DT)
GRID_Z = numpy.interp(GRID_TIME, 
                      universe_age, 
                      universe_redshifts)





class galaxy(object):

    def __init__(self):
        self.tau = 0.        # tau for initail exp-declining SFH
        self.dt = 0.         # discretized time unit (Gyr)
        self.sfr = numpy.array([])
        self.stellar_mass = numpy.array([])
        self.UV_color = numpy.array([])
        self.VJ_color = numpy.array([])
        self.sfqu_flag = numpy.array([])
        self.cpu_time = 0.    # time it took to calculate this galaxy (seconds)
        self.GRID_TIME = GRID_TIME
        self.GRID_Z = GRID_Z

    def append_stellar_mass(self, m):
        self.stellar_mass = numpy.append(self.stellar_mass, m)

    def append_UV_color(self, uv):
        self.UV_color = numpy.append(self.UV_color, uv)

    def append_VJ_color(self, vj):
        self.VJ_color = numpy.append(self.VJ_color, vj)

    def append_sfqu_flag(self, sfqu):
        self.sfqu_flag = numpy.append(self.sfqu_flag, sfqu)

    def get_sfr(self, tage):
        '''
        Returns the SFR at the requested age (in Gyr)
        '''
        return numpy.interp(tage, self.GRID_TIME, self.sfr)

    def get_ssfr(self, tage):
        '''
        Returns the SSFR at the requested age (in Gyr)
        '''
        return numpy.interp(tage, self.GRID_TIME, self.sfr/self.stellar_mass)

    def get_stellar_mass(self, tage):
        '''
        Returns the stellar mass at the requested age (in Gyr)
        '''
        return numpy.interp(tage, self.GRID_TIME, self.stellar_mass)

    def get_UV_color(self, tage):
        '''
        Returns the rest-frame (U-V) at the requested age (in Gyr)
        '''
        return numpy.interp(tage, self.GRID_TIME, self.UV_color)

    def get_VJ_color(self, tage):
        '''
        Returns the rest-frame (V-J) at the requested age (in Gyr)
        '''
        return numpy.interp(tage, self.GRID_TIME, self.VJ_color)

    def get_sfqu_flag(self, tage):
        '''
        Returns the SF/QU flag at the requested age (in Gyr)
        '''
        uv = self.get_UV_color(tage)
        vj = self.get_VJ_color(tage)
        z = numpy.interp(tage, self.GRID_TIME, self.GRID_Z)
        return mypy.uvj_select(uv, vj, [z])[0]



def read_galaxy_file(filename):

    g = galaxy()

    fopen = open(filename, 'r')
    line = '#'
    while line[0] == '#':
        line = fopen.readline()
        split = line.split()
        if split[1] == 'tau':
            g.tau = float(split[3])
        if split[1] == 'tform_best':
            g.tform_best = float(split[3])
        if split[1] == 'cpu_time':
            g.cpu_time = float(split[3])
    fopen.close()

    data = numpy.loadtxt(filename)
    g.sfr = data[:,0]
    g.stellar_mass = data[:,1]
    g.UV_color = data[:,2]
    g.VJ_color = data[:,3]

    del data
    return g






########################################################
###  Reading in simulated galaxies from ascii files  ###
########################################################

if True:

    zstart = 1.05
    zend = 0.8

    tstart = cosmo.age(zstart).value
    tend = cosmo.age(zend).value


    dir0 = '../data/simulations/single_exp'
    dir0 = '/Users/atomczak/tmp_transfer/single_exp'


    sims_data = []

    sims_data.append([glob.glob('%s/0.2000/galaxy*txt' % dir0), 0.20, '0.2'])
    sims_data.append([glob.glob('%s/0.3000/galaxy*txt' % dir0), 0.30, '0.3'])
    sims_data.append([glob.glob('%s/0.4000/galaxy*txt' % dir0), 0.40, '0.4'])
    sims_data.append([glob.glob('%s/0.5000/galaxy*txt' % dir0), 0.50, '0.5'])
    sims_data.append([glob.glob('%s/0.7500/galaxy*txt' % dir0), 0.75, '0.75'])
    sims_data.append([glob.glob('%s/1.0000/galaxy*txt' % dir0), 1.00, '1'])
    sims_data.append([glob.glob('%s/2.0000/galaxy*txt' % dir0), 2.00, '2'])
    sims_data.append([glob.glob('%s/4.0000/galaxy*txt' % dir0), 4.00, '4'])
    sims_data.append([glob.glob('%s/8.0000/galaxy*txt' % dir0), 8.00, '8'])

    filenames_all = [sd[0] for sd in sims_data]
    taus     = numpy.array([sd[1] for sd in sims_data])
    taus_str = numpy.array([sd[2] for sd in sims_data])

    my_galaxies = []

    for i_tau in range(len(taus)):

        print('reading tau = %.4f Gyr files' % taus[i_tau])

        my_galaxies.append([])
        N = len(filenames_all[i_tau])
        ###  Initializing progress bar
        widgets = ['  Running: ', progressbar.Percentage(), 
                   ' ', progressbar.Bar(marker=progressbar.RotatingMarker()), 
                   ' ', progressbar.ETA(), 
                   ' ', progressbar.FileTransferSpeed()]
        pbar = progressbar.ProgressBar(widgets=widgets, maxval=N).start()

        for i_filename in range(N):

            pbar.update(i_filename)

            f = filenames_all[i_tau][i_filename]
            my_galaxies[-1].append(read_galaxy_file(f))
        
        print('\n')






#############################################################
###  Plotting median SSFR in M* bins between zstart/zend  ###
###  Single panel
#############################################################

if False:

    lmass_start = []
    lsfr_start  = []
    lmass_end   = []
    lsfr_end    = []

    for i_gal in range(len(my_galaxies)):

        g = my_galaxies[i_gal]

        if g.get_sfqu_flag(tstart) == 1:
            lmass_start.append(numpy.log10(g.get_stellar_mass(tstart)))
            lsfr_start.append(numpy.log10(g.get_sfr(tstart)))

        if g.get_sfqu_flag(tend) == 1:
            lmass_end.append(numpy.log10(g.get_stellar_mass(tend)))
            lsfr_end.append(numpy.log10(g.get_sfr(tend)))

    lmass_start = numpy.array(lmass_start)
    lsfr_start  = numpy.array(lsfr_start)
    lmass_end   = numpy.array(lmass_end)
    lsfr_end    = numpy.array(lsfr_end)

    lssfr_start = lsfr_start - lmass_start
    lssfr_end = lsfr_end - lmass_end



    ###  binning the data
    lmassbins = numpy.array([9.4, 10.1, 10.8])
    lmassbin_labels = [' %4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$) < %4.1f' % (lmassbins[0], lmassbins[1]),
                       '%4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$) < %4.1f' % (lmassbins[1], lmassbins[2]),
                       '%4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$)' % lmassbins[2]]

    digi_lmass_start = numpy.digitize(lmass_start, lmassbins)
    digi_lmass_end = numpy.digitize(lmass_end, lmassbins)




    ###  setting up figure
    fig_ssfr = pyplot.figure(figsize=(7.36, 7.17))

    sp1_ssfr = fig_ssfr.add_subplot(111)
    upper_axis = sp1_ssfr.twiny()

    sp1_ssfr.grid(color='#666666')
    sp1_ssfr.minorticks_on()
    sp1_ssfr.set_xlabel('age of Universe [Gyr]', size=22)
    sp1_ssfr.set_ylabel('log( SSFR$_{UV+IR}$ / [yr$^{-1}$] )', size=22)
    upper_axis.set_xlabel('redshift', size=20)

    dt_baseline = abs(tend - tstart)
    t_axis_1 = tstart-dt_baseline / 3.
    t_axis_2 = tend+dt_baseline / 3.
    z_axis_1 = numpy.interp(t_axis_1, universe_age, universe_redshifts)
    z_axis_2 = numpy.interp(t_axis_2, universe_age, universe_redshifts)

    sp1_ssfr.axis([t_axis_1, t_axis_2, -10.12, -8.58])
    upper_axis.set_xlim(z_axis_1, z_axis_2)

    fig_ssfr.subplots_adjust(wspace=0, left=0.17, right=0.98)
    fig_ssfr.subplots_adjust(wspace=0, left=0.17, right=0.98, top=0.92)


    colors = ['b', 'yellow', 'orange']

    dt02 = 0.02 * dt_baseline
    offies = [-dt02, 0, dt02]

    for i_digi_lmass in [1, 2, 3]:

        inds_start = numpy.where(digi_lmass_start == i_digi_lmass)[0]
        inds_end = numpy.where(digi_lmass_end == i_digi_lmass)[0]

        print(len(inds_start), len(inds_end))

        p_16_50_84_start = numpy.percentile(lssfr_start[inds_start], [16, 50, 84])
        p_16_50_84_end = numpy.percentile(lssfr_end[inds_end], [16, 50, 84])

        p16_start, p50_start, p84_start = p_16_50_84_start
        p16_end, p50_end, p84_end = p_16_50_84_end

        elo_start = p50_start - p16_start
        ehi_start = p84_start - p50_start
        elo_end = p50_end - p16_end
        ehi_end = p84_end - p50_end

        offie = offies[i_digi_lmass-1]

        sp1_ssfr.plot([tstart+offie, tend+offie], 
                      [p50_start, p50_end],
                      color=colors[i_digi_lmass-1], 
                      lw=2, 
                      zorder=1,
                      label=lmassbin_labels[i_digi_lmass-1],
                      path_effects=[patheffects.Stroke(linewidth=6, 
                                                       foreground='k'), 
                                    patheffects.Normal()])

        sp1_ssfr.plot([tstart+offie, tend+offie], 
                      [p50_start, p50_end],
                      ls='',
                      marker='o',
                      ms=9,
                      mew=2,
                      mfc=colors[i_digi_lmass-1], 
                      zorder=2)

    sp1_ssfr.legend(loc=3, numpoints=1, fontsize=16)


    ###  Case 1: single exponential
    if True:
        title = 'Model galaxies\n'
        title += 'SFH = single exp.\n'
        title += r'$\tau$ = 2 Gyr'

    ###  Case 2a: double exponential, t_trans = 0.5 Gyr
    if False:
        title = 'Model galaxies\n'
        title += 'SFH = double exp.\n'
        title += r'$\tau_1$ = 2 Gyr'
        title += '\n'
        title += r'$\tau_2$ = 0.5 Gyr'
        title += '\n'
        title += r'$t_{trans}$ = 0.5 Gyr'


    t = sp1_ssfr.text(0.97, 
                      0.97, 
                      title, 
                      fontsize=16,
                      transform=sp1_ssfr.transAxes,
                      horizontalalignment='right',
                      verticalalignment='top')












#############################################################
###  Plotting median SSFR in M* bins between zstart/zend  ###
###  Multi panel
#############################################################

if True:

    fig = pyplot.figure(figsize=(14.7, 10.3))

    sp3 = fig.add_subplot(233)
    sp2 = fig.add_subplot(232)
    sp1 = fig.add_subplot(231)
    sp6 = fig.add_subplot(236)
    sp5 = fig.add_subplot(235)
    sp4 = fig.add_subplot(234)

    sps_top = [sp1, sp2, sp3]
    sps_bottom = [sp4, sp5, sp6]

    fig.subplots_adjust(wspace=0, left=0.09, bottom=0.09, top=0.99, right=0.99)


    for sp in sps_top + sps_bottom:
        sp.grid()
        sp.minorticks_on()

    for sp in sps_top:
        sp.set_xlabel('log( 1 + $\delta_{gal}$ )')
        sp.set_ylabel('log( SSFR / [yr$^{-1}$] )')
        sp.axis([0.05, 1.44, -10.1, -8.7])

    for sp in sps_bottom:
        sp.set_xlabel('redshift')
        sp.set_ylabel('log( SSFR / [yr$^{-1}$] )')
        sp.axis([1.17, 0.69, -10.1, -8.7])

        sp.xaxis.set_ticks(numpy.arange(0.7, 1.15, 0.1))






    ###  reading in ORELSE data
    colors = ['b', 'yellow', 'orange']
    lmassbins = numpy.array([9.4, 10.1, 10.8])
    lmassbin_labels = [' %4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$) < %4.1f' % (lmassbins[0], lmassbins[1]),
                       '%4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$) < %4.1f' % (lmassbins[1], lmassbins[2]),
                       '%4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$)' % lmassbins[2]]
    lmassbin_labels = [u'log($\mathit{M}$): %.1f\u2212%.1f'  % (lmassbins[0], lmassbins[1]),
                       u'log($\mathit{M}$): %.1f\u2212%.1f'  % (lmassbins[1], lmassbins[2]),
                       u'log($\mathit{M}$): >%.1f'  % lmassbins[2]]

    overdens_median_lowz_sf = numpy.loadtxt('../data/table_SSFR_density_Mstar_lowz___overdens_median_sf.dat')
    lssfr_median_lowz_sf = numpy.loadtxt('../data/table_SSFR_density_Mstar_lowz___lssfr_median_sf.dat')
    lsfr_median_lowz_sf_err = numpy.loadtxt('../data/table_SSFR_density_Mstar_lowz___lsfr_median_sf_err.dat')

    overdens_median_highz_sf = numpy.loadtxt('../data/table_SSFR_density_Mstar_highz___overdens_median_sf.dat')
    lssfr_median_highz_sf = numpy.loadtxt('../data/table_SSFR_density_Mstar_highz___lssfr_median_sf.dat')
    lsfr_median_highz_sf_err = numpy.loadtxt('../data/table_SSFR_density_Mstar_highz___lsfr_median_sf_err.dat')


    ###  plotting ORELSE data
    for i_sp in range(len(sps_top)):

        sp = sps_top[i_sp]

        p1 = sp.errorbar(overdens_median_lowz_sf[i_sp],
                         lssfr_median_lowz_sf[i_sp],
                         yerr=lsfr_median_lowz_sf_err[i_sp],
                         dashes=[5,2], lw=2, marker='v', mew=1.8, ms=9, elinewidth=2, zorder=2,
                         mfc=colors[i_sp], color='k', mec='k', ecolor='k')

        p2 = sp.errorbar(overdens_median_highz_sf[i_sp],
                         lssfr_median_highz_sf[i_sp],
                         yerr=lsfr_median_highz_sf_err[i_sp],
                         ls='-', lw=2, marker='^', mew=1.8, ms=9, elinewidth=2, zorder=2,
                         mfc=colors[i_sp], color='k', mec='k', ecolor='k')
        t = sp.text(0.05, 0.06, lmassbin_labels[i_sp], transform=sp.transAxes,
                    horizontalalignment='left', verticalalignment='bottom')

    for i_sp in range(len(sps_bottom)):
        sp = sps_bottom[i_sp]
        t = sp.text(0.05, 0.06, lmassbin_labels[i_sp], transform=sp.transAxes,
                    horizontalalignment='left', verticalalignment='bottom')



    ###  adding legend data
    sp3.errorbar(0, 0, yerr=[0],
                 ls='-', lw=2, marker='^', mew=1.8, ms=9, elinewidth=2, zorder=2,
                 mfc='w', color='k', mec='k', ecolor='k', 
                 label='0.9 < z < 1.3')
    sp3.errorbar(0, 0, yerr=[0],
                 dashes=[5,2], lw=2, marker='v', mew=1.8, ms=9, elinewidth=2, zorder=2,
                 mfc='w', color='k', mec='k', ecolor='k', 
                 label='0.6 < z < 0.9')
    leg = sp3.legend(loc=2, fontsize=16)





    ###  plotting simulation data

    plot_lines = []
    lssfrs = []

    for i_tau in range(len(taus)):

        plot_lines.append([])
        lssfrs.append([])
        galaxies = my_galaxies[i_tau]

        ###  getting quantities at zstart=1.05 and zend=0.8
        galaxies_lmass_start = numpy.array([numpy.log10(g.get_stellar_mass(tstart)) for g in galaxies])
        galaxies_lmass_end = numpy.array([numpy.log10(g.get_stellar_mass(tend)) for g in galaxies])

        galaxies_lssfr_start = numpy.array([numpy.log10(g.get_ssfr(tstart)) for g in galaxies])
        galaxies_lssfr_end = numpy.array([numpy.log10(g.get_ssfr(tend)) for g in galaxies])

        galaxies_sfqu_start = numpy.array([g.get_sfqu_flag(tstart) for g in galaxies])
        galaxies_sfqu_end = numpy.array([g.get_sfqu_flag(tend) for g in galaxies])


        digi_lmass_start = numpy.digitize(galaxies_lmass_start, lmassbins)
        digi_lmass_end = numpy.digitize(galaxies_lmass_end, lmassbins)

        for i_sp in range(len(sps_bottom)):

            sp = sps_bottom[i_sp]

            ###  identifying galaxies in appropriate M* bin at zstart and zend
            inds_lmass_start = numpy.where((digi_lmass_start == i_sp+1) & (galaxies_sfqu_start == 1))[0]
            inds_lmass_end = numpy.where((digi_lmass_end == i_sp+1) & (galaxies_sfqu_end == 1))[0]

            lssfr_median_highz_sf_sim = numpy.median(galaxies_lssfr_start[inds_lmass_start])
            lssfr_median_lowz_sf_sim = numpy.median(galaxies_lssfr_end[inds_lmass_end])

            pl = sp.plot([zstart, zend],
                         [lssfr_median_highz_sf_sim, lssfr_median_lowz_sf_sim],
                         lw=3,
                         color=pyplot.cm.rainbow_r(i_tau * 1. / (len(taus) - 1.)),
                         label=r'$\tau$ = %s Gyr' % taus_str[i_tau])
                         #label=r'$\tau$ = %.2f Gyr' % taus[i_tau])

            plot_lines[-1].append(pl[0])
            lssfrs[-1].append([lssfr_median_highz_sf_sim, lssfr_median_lowz_sf_sim])

    lssfrs = numpy.array(lssfrs)


    ###  shifting curves to mean lssfr at high-z
    for i_sp in range(len(sps_bottom)):

        lssfr_zeropoint = numpy.average(lssfrs[:, i_sp, 0])
        offsets = lssfrs[:, i_sp, 0] - lssfr_zeropoint

        for i_tau in range(len(taus)):

            plot_lines[i_tau][i_sp].set_ydata(lssfrs[i_tau, i_sp] - offsets[i_tau])



    leg6 = sp6.legend(loc=2, fontsize=14, ncol=2)

    pyplot.savefig('../figures/simulation_results2.pdf')
    pyplot.close()




#################################################
###  Creating diagnostic plots for simulated
###  galaxies at zstart and zend
#################################################

if False:


    pdf = PdfPages('../figures/simulation_diagnostic_plots.pdf')

    fig = pyplot.figure(figsize=(14.7, 7.2))

    sp2 = fig.add_subplot(122)
    sp1 = fig.add_subplot(121)

    fig.subplots_adjust(wspace=0, left=0.1)


    colors = ['b', 'yellow', 'orange']
    lmassbins = numpy.array([9.4, 10.1, 10.8])
    lmassbin_labels = [' %4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$) < %4.1f' % (lmassbins[0], lmassbins[1]),
                       '%4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$) < %4.1f' % (lmassbins[1], lmassbins[2]),
                       '%4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$)' % lmassbins[2]]


    for i_tau in range(len(taus))[::-1]:

        sp1.minorticks_on()
        sp2.minorticks_on()

        sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
        sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
        sp1.set_ylabel('log( SFR / [M$_{\odot}$ / yr$^{-1}$] )')
        sp2.set_ylabel('log( SFR / [M$_{\odot}$ / yr$^{-1}$] )')

        sp1.axis([8.3, 11.7, -0.6, 2.2])
        sp2.axis([8.3, 11.7, -0.6, 2.2])


        galaxies = my_galaxies[i_tau]

        ###  getting quantities at zstart=1.05 and zend=0.8
        galaxies_lmass_start = numpy.array([numpy.log10(g.get_stellar_mass(tstart)) for g in galaxies])
        galaxies_lmass_end = numpy.array([numpy.log10(g.get_stellar_mass(tend)) for g in galaxies])

        galaxies_lsfr_start = numpy.array([numpy.log10(g.get_sfr(tstart)) for g in galaxies])
        galaxies_lsfr_end = numpy.array([numpy.log10(g.get_sfr(tend)) for g in galaxies])

        galaxies_sfqu_start = numpy.array([g.get_sfqu_flag(tstart) for g in galaxies])
        galaxies_sfqu_end = numpy.array([g.get_sfqu_flag(tend) for g in galaxies])


        ###  counting galaxies in massbins
        digi_lmass_start = numpy.digitize(galaxies_lmass_start, lmassbins)
        digi_lmass_end = numpy.digitize(galaxies_lmass_end, lmassbins)

        inds_sf_start = numpy.where(galaxies_sfqu_start == 1)[0]
        inds_sf_end = numpy.where(galaxies_sfqu_end == 1)[0]

        inds_qu_start = numpy.where(galaxies_sfqu_start == 0)[0]
        inds_qu_end = numpy.where(galaxies_sfqu_end == 0)[0]

        ntot_start = numpy.bincount(digi_lmass_start)
        ntot_end = numpy.bincount(digi_lmass_end)
        nqu_start = ntot_start * 0.
        nqu_end = ntot_end * 0.


        ###  plotting SF and QU galaxies
        sp1.plot(galaxies_lmass_start[inds_sf_start], galaxies_lsfr_start[inds_sf_start], 'ko', ms=1, mec='b')
        sp1.plot(galaxies_lmass_start[inds_qu_start], galaxies_lsfr_start[inds_qu_start], 'ko', ms=1, mec='r')
        sp2.plot(galaxies_lmass_end[inds_sf_end], galaxies_lsfr_end[inds_sf_end], 'ko', ms=1, mec='b')
        sp2.plot(galaxies_lmass_end[inds_qu_end], galaxies_lsfr_end[inds_qu_end], 'ko', ms=1, mec='r')

        label1 = 'z = 1.05\n'
        label1 += r'$\tau$ = %.2f' % taus[i_tau]
        label1 += '\nUVJ Star-forming'
        t1 = sp1.text(0.03, 0.95, label1, transform=sp1.transAxes, fontsize=20, 
                      horizontalalignment='left', verticalalignment='top')

        label2 = 'z = 0.8\n'
        label2 += r'$\tau$ = %.2f' % taus[i_tau]
        label2 += '\nUVJ Star-forming'
        t2 = sp2.text(0.03, 0.95, label2, transform=sp2.transAxes, fontsize=20, 
                      horizontalalignment='left', verticalalignment='top')



        ###  plotting numbers of galaxies in M* bins

        if len(inds_qu_start) > 0:
            digi_lmass_qu_start = numpy.digitize(galaxies_lmass_start[inds_qu_start], lmassbins)
            nqu_start = numpy.bincount(digi_lmass_qu_start)
        if len(inds_qu_end) > 0:
            digi_lmass_qu_end = numpy.digitize(galaxies_lmass_end[inds_qu_end], lmassbins)
            nqu_end = numpy.bincount(digi_lmass_qu_end)


        for i_lmass in range(len(lmassbins)):
            sp1.axvline(lmassbins[i_lmass], color='gray', lw=1)
            sp2.axvline(lmassbins[i_lmass], color='gray', lw=1)

            ngals1 = 'QUI: %i' % nqu_start[i_lmass+1]
            ngals1 += '\nTOT: %i' % ntot_start[i_lmass+1]
            ta = sp1.text(lmassbins[i_lmass], 1.95, ngals1, fontsize=13)

            ngals2 = 'QUI: %i' % nqu_end[i_lmass+1]
            ngals2 += '\nTOT: %i' % ntot_end[i_lmass+1]
            tb = sp2.text(lmassbins[i_lmass], 1.95, ngals2, fontsize=13)



        pdf.savefig(fig)
        sp1.clear()
        sp2.clear()

    pdf.close()
    pyplot.close()









#################################################
###  Creating schematic figure that illustrates
###  how the simulation works
#################################################

if True:

    fig = pyplot.figure(figsize=(12.3, 10.3))

    sp0 = pyplot.subplot2grid((7, 8), (0, 1), rowspan=3, colspan=6)
    sp2 = pyplot.subplot2grid((7, 8), (3, 4), rowspan=4, colspan=4)
    sp1 = pyplot.subplot2grid((7, 8), (3, 0), rowspan=4, colspan=4)

    fig.subplots_adjust(wspace=0.13, hspace=1., left=0.09, bottom=0.08, top=0.94)

    colors = ['b', 'yellow', 'orange']
    lmassbins = numpy.array([9.4, 10.1, 10.8])
    lmassbin_labels = [' %4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$) < %4.1f' % (lmassbins[0], lmassbins[1]),
                       '%4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$) < %4.1f' % (lmassbins[1], lmassbins[2]),
                       '%4.1f < log($\mathit{M}_*$/$\mathit{M}_{\odot}$)' % lmassbins[2]]






    sp0.minorticks_on()
    sp1.minorticks_on()
    sp2.minorticks_on()

    sp0.set_xlabel('age of Universe [Gyr]')
    sp0.set_ylabel('log( SFR / [$\mathit{M}_{\odot}$ / yr$^{-1}$] )')

    sp1.set_xlabel('log( $\mathit{M}_*$ / $\mathit{M}_{\odot}$ )')
    sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp1.set_ylabel('log( SFR / [$\mathit{M}_{\odot}$ / yr$^{-1}$] )')
    sp2.set_ylabel('log( SFR / [$\mathit{M}_{\odot}$ / yr$^{-1}$] )')

    sp1.axis([8.3, 11.7, -0.6, 2.2])
    sp2.axis([8.3, 11.7, -0.6, 2.2])
    sp2.yaxis.set_ticklabels([])






    ###  plotting tau=0.75 simulation
    i_tau = taus.tolist().index(0.75)
    galaxies = my_galaxies[i_tau]

    ###  getting quantities at zstart=1.05 and zend=0.8
    galaxies_lmass_start = numpy.array([numpy.log10(g.get_stellar_mass(tstart)) for g in galaxies])
    galaxies_lmass_end = numpy.array([numpy.log10(g.get_stellar_mass(tend)) for g in galaxies])

    galaxies_lsfr_start = numpy.array([numpy.log10(g.get_sfr(tstart)) for g in galaxies])
    galaxies_lsfr_end = numpy.array([numpy.log10(g.get_sfr(tend)) for g in galaxies])

    galaxies_sfqu_start = numpy.array([g.get_sfqu_flag(tstart) for g in galaxies])
    galaxies_sfqu_end = numpy.array([g.get_sfqu_flag(tend) for g in galaxies])


    ###  counting galaxies in massbins
    digi_lmass_start = numpy.digitize(galaxies_lmass_start, lmassbins)
    digi_lmass_end = numpy.digitize(galaxies_lmass_end, lmassbins)

    inds_sf_start = numpy.where(galaxies_sfqu_start == 1)[0]
    inds_sf_end = numpy.where(galaxies_sfqu_end == 1)[0]

    inds_qu_start = numpy.where(galaxies_sfqu_start == 0)[0]
    inds_qu_end = numpy.where(galaxies_sfqu_end == 0)[0]

    ntot_start = numpy.bincount(digi_lmass_start)
    ntot_end = numpy.bincount(digi_lmass_end)
    nqu_start = ntot_start * 0.
    nqu_end = ntot_end * 0.


    ###  plotting SF and QU galaxies
    sp1.plot(galaxies_lmass_start[inds_sf_start], galaxies_lsfr_start[inds_sf_start], 'ko', ms=1, mec='b')
    #sp1.plot(galaxies_lmass_start[inds_qu_start], galaxies_lsfr_start[inds_qu_start], 'ko', ms=1, mec='r')
    sp2.plot(galaxies_lmass_end[inds_sf_end], galaxies_lsfr_end[inds_sf_end], 'ko', ms=1, mec='b')
    sp2.plot(galaxies_lmass_end[inds_qu_end], galaxies_lsfr_end[inds_qu_end], 'ko', ms=1, mec='r')

    label1 = 'z = %.2f\n' % zstart
    label1 += r'$\tau$ = %.2f Gyr' % taus[i_tau]
    t1 = sp1.text(0.04, 0.95, label1, transform=sp1.transAxes, fontsize=17, 
                  horizontalalignment='left', verticalalignment='top')

    label2 = 'z = %.2f\n' % zend
    label2 += r'$\tau$ = %.2f Gyr' % taus[i_tau]
    t2 = sp2.text(0.04, 0.95, label2, transform=sp2.transAxes, fontsize=17, 
                  horizontalalignment='left', verticalalignment='top')



    ###  plotting numbers of galaxies in M* bins
    if len(inds_qu_start) > 0:
        digi_lmass_qu_start = numpy.digitize(galaxies_lmass_start[inds_qu_start], lmassbins)
        nqu_start = numpy.bincount(digi_lmass_qu_start)
    if len(inds_qu_end) > 0:
        digi_lmass_qu_end = numpy.digitize(galaxies_lmass_end[inds_qu_end], lmassbins)
        nqu_end = numpy.bincount(digi_lmass_qu_end)


    for i_lmass in range(len(lmassbins)):
        sp1.axvline(lmassbins[i_lmass], color='gray', lw=1)
        sp2.axvline(lmassbins[i_lmass], color='gray', lw=1)

        ngals1 = 'QUI: %i' % nqu_start[i_lmass+1]
        ngals1 += '\nTOT: %i' % ntot_start[i_lmass+1]
        #ta = sp1.text(lmassbins[i_lmass], 1.95, ngals1, fontsize=13)

        ngals2 = 'QUI: %i' % nqu_end[i_lmass+1]
        ngals2 += '\nTOT: %i' % ntot_end[i_lmass+1]
        #tb = sp2.text(lmassbins[i_lmass], 1.95, ngals2, fontsize=13)




    ###  plotting upper panel with SFHs

    sp0.axvline(tstart, color='k', lw=2, dashes=[5, 2])
    sp0.axvline(tend, color='k', lw=2, dashes=[5, 2])

    x_offie = 0.04
    t01 = sp0.text(tstart+x_offie, 1.41, 'z = %.2f' % zstart, 
                   fontsize=13, verticalalignment='top')
    t02 = sp0.text(tend+x_offie, 1.41, 'z = %.2f' % zend, 
                   fontsize=13, verticalalignment='top')





    ###  finding 1 galaxy in each tau realization with SFR(tstart) ~ 10 Msol/yr
    plot_gals = []
    for i_tau in range(len(taus)):
        if taus[i_tau] < 0.5: continue
        for i_gal in range(len(my_galaxies[i_tau])):
            g = my_galaxies[i_tau][i_gal]
            sfr_tstart = g.get_sfr(tstart)
            lmass_tstart = numpy.log10(g.get_stellar_mass(tstart))
            if 4.75 < sfr_tstart < 5.25 and 9.775 < lmass_tstart < 9.825:
                plot_gals.append(g)
                sp0.plot(GRID_TIME, numpy.log10(g.sfr+1.e-10), 
                         color=pyplot.cm.rainbow_r(i_tau * 1. / (len(taus) - 1.)))



                lsfr_max = numpy.log10(g.sfr.max())
                tform = GRID_TIME[numpy.where(g.sfr == g.sfr.max())[0][0]]

                t_tau = sp0.text(tform, lsfr_max, taus_str[i_tau], fontsize=15,
                                 color=pyplot.cm.rainbow_r(i_tau * 1. / (len(taus) - 1.)),
                                 horizontalalignment='right', verticalalignment='bottom')

                break


    sp0.axis([3.6, 7.7, -0.3, 1.5])


    sp0.xaxis.set_label_position('top')
    sp0.xaxis.set_ticks_position('top')
    sp0.xaxis.set_ticks_position('both')
    fig.subplots_adjust(hspace=0.4)


    pyplot.savefig('../figures/simulation_schematic2.pdf')
    pyplot.close()


    #top_axis = sp0.twiny()
    #top_axis.set_xlabel('redshift')

print '\nDONE!!!\n'









#########################################
###  Plotting quiescent fraction vs. M*
###  for simulated galaxies at z=0.8
#########################################

if True:

    dm = 0.7
    lmassbins = numpy.array([9.4, 10.1, 10.8])
    lmassbars = numpy.array([9.75, 10.45, 11.15])

    ###  numbers of SF and TOT galaxies in each M* bin for each tau realization of simulation at tstart and tend
    N_SF_tstart  = numpy.zeros((len(sims_data), len(lmassbins)), dtype=float)
    N_SF_tend    = numpy.zeros((len(sims_data), len(lmassbins)), dtype=float)
    N_TOT_tstart = numpy.zeros((len(sims_data), len(lmassbins)), dtype=float)
    N_TOT_tend   = numpy.zeros((len(sims_data), len(lmassbins)), dtype=float)
    lmasses_TOT_tstart = numpy.zeros((len(sims_data), len(lmassbins)), dtype=float)
    lmasses_TOT_tend   = numpy.zeros((len(sims_data), len(lmassbins)), dtype=float)

    colors = pyplot.cm.rainbow_r(numpy.linspace(0., 1., len(sims_data)))


    for i_tau in range(len(sims_data)):

        tau     = sims_data[i_tau][1]
        tau_str = sims_data[i_tau][2]

        lmass_gals_tstart = numpy.array([numpy.log10(g.get_stellar_mass(tstart)) for g in my_galaxies[i_tau]])
        sfqu_gals_tstart  = numpy.array([g.get_sfqu_flag(tstart) for g in my_galaxies[i_tau]])

        lmass_gals_tend = numpy.array([numpy.log10(g.get_stellar_mass(tend)) for g in my_galaxies[i_tau]])
        sfqu_gals_tend  = numpy.array([g.get_sfqu_flag(tend) for g in my_galaxies[i_tau]])


        ###  digitizing galaxies into M* bins
        digi_lmass_tstart = numpy.digitize(lmass_gals_tstart, lmassbins)
        digi_lmass_tend = numpy.digitize(lmass_gals_tend, lmassbins)

        for i_lmass in range(len(lmassbins)):

            inds_lmass_tstart = numpy.where(digi_lmass_tstart == i_lmass+1)[0]
            inds_lmass_tend = numpy.where(digi_lmass_tend == i_lmass+1)[0]

            ntot_tstart = len(inds_lmass_tstart)
            ntot_tend   = len(inds_lmass_tend)
            nsf_tstart  = numpy.count_nonzero(sfqu_gals_tstart[inds_lmass_tstart])
            nsf_tend    = numpy.count_nonzero(sfqu_gals_tend[inds_lmass_tend])

            N_SF_tstart[i_tau][i_lmass]  = nsf_tstart
            N_SF_tend[i_tau][i_lmass]  = nsf_tend
            N_TOT_tstart[i_tau][i_lmass] = ntot_tstart
            N_TOT_tend[i_tau][i_lmass] = ntot_tend
            lmasses_TOT_tstart[i_tau][i_lmass] = numpy.median(lmass_gals_tstart[inds_lmass_tstart])
            lmasses_TOT_tend[i_tau][i_lmass] = numpy.median(lmass_gals_tend[inds_lmass_tend])

    N_QU_tstart = N_TOT_tstart - N_SF_tstart
    N_QU_tend   = N_TOT_tend - N_SF_tend
    f_newly_quenched_tstart = N_QU_tstart / N_TOT_tstart
    f_newly_quenched_tend   = N_QU_tend / N_TOT_tend




    ###  plotting "newly quenched" fractions 
    ###  (faction of galaxies which began as SF at z=1.05 and quenched at z=0.8)

    fig = pyplot.figure()
    sp = fig.add_subplot(111)

    sp.grid()
    sp.minorticks_on()
    sp.set_xlabel('$log( M_* / M_{\odot} )$')
    #sp.set_ylabel('Quiescent Fraction')
    sp.set_ylabel("''Newly-Quenched'' Fraction")
    sp.axis([9.3, 11.7, -0.1, 1.1])


    ###  all taus in one panel
    for i_tau in range(len(sims_data)):

        sp.errorbar(lmasses_TOT_tend[i_tau], f_newly_quenched_tend[i_tau], 
                    color=colors[i_tau], lw=3, marker='o', 
                    ms=8, ecolor='k', elinewidth=1.5, mew=2.,
                    label=taus_str[i_tau])

    sp.legend(loc=1, fontsize=14, title=r'$\tau$ / Gyr')


    fig.savefig('../figures/simulation_newly_quenched_fraction.pdf')







    ###  grabbing quiescent fractions from ZFOURGE at z=1.05 and 0.8
    fqu_zfourge_tstart = fquiescent_leja2015(1.05, lmasses_TOT_tstart)
    fqu_zfourge_tend   = fquiescent_leja2015(0.8, lmasses_TOT_tend)


    ###  Using these quiescent fractions, we calculate the number of 
    ###  quiescent galaxies at z=1.05 that are not being considered in
    ###  the simulation. Assuming this already existing population
    ###  of quiescent galaxies simply carry on to z=0.8 (i.e. without
    ###  rejuvenation) we can estimate the quiescent fraction at z=0.8
    ###  for each tau realization of the simulation.

    N_QU_existing_tstart = N_SF_tstart * fqu_zfourge_tstart / (1 - fqu_zfourge_tstart)



    ###  Now that we've estimated the number of "missing" quiescent
    ###  galaxies at z=1.05, we can carry them over to z=0.8 to infer
    ###  the quiescent fraction for each tau realization

    fqu_tstart = (N_QU_tstart + N_QU_existing_tstart) / (N_TOT_tstart + N_QU_existing_tstart)
    fqu_tend   = (N_QU_tend + N_QU_existing_tstart) / (N_TOT_tend + N_QU_existing_tstart)











    ###  plotting 4-panel figure, to illustrate process

    fig = pyplot.figure(figsize=(15.4, 13.))

    sp2 = fig.add_subplot(222)
    sp1 = fig.add_subplot(221)
    sp4 = fig.add_subplot(224)
    sp3 = fig.add_subplot(223)

    sp1.minorticks_on()
    sp2.minorticks_on()
    sp3.minorticks_on()
    sp4.minorticks_on()

    sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp3.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp4.set_xlabel('log( M$_*$ / M$_{\odot}$ )')

    sp1.set_ylabel('Quiescent Fraction')
    sp2.set_ylabel('Quiescent Fraction')
    sp3.set_ylabel('Quiescent Fraction')
    sp4.set_ylabel('Quiescent Fraction')

    sp1.set_title('z = 1.05')
    sp2.set_title('z = 0.8')

    sp1.axis([9.3, 11.6, -0.1, 1.19])
    sp2.axis([9.3, 11.6, -0.1, 1.19])
    sp3.axis([9.3, 11.6, -0.1, 1.19])
    sp4.axis([9.3, 11.6, -0.1, 1.19])

    sp1.grid()
    sp2.grid()
    sp3.grid()
    sp4.grid()

    fig.subplots_adjust(hspace=0, left=0.07, bottom=0.07)




    text2 = '"Newly-Quenched" from\nz=1.05 SF galaxies'
    t2 = sp2.text(0.03, 0.97, text2, transform=sp2.transAxes,
                  horizontalalignment='left', verticalalignment='top')

    text3 = 'with QU galaxies added\nto match ZFOURGE at z=1.05'
    t3 = sp3.text(0.03, 0.97, text3, transform=sp3.transAxes,
                  horizontalalignment='left', verticalalignment='top')

    text4 = 'Final predictions'
    t4 = sp4.text(0.03, 0.97, text4, transform=sp4.transAxes,
                  horizontalalignment='left', verticalalignment='top')


    linestyles = ['-', '--']

    for i_tau in range(len(taus)):

        ###  plotting initial quiescent fractions
        sp1.errorbar(lmasses_TOT_tstart[i_tau], f_newly_quenched_tstart[i_tau],
                     color=colors[i_tau], lw=4, marker='o', ls=linestyles[i_tau%2], 
                     ms=8, ecolor='k', elinewidth=1.5, mew=2.,
                     label=taus_str[i_tau])


        ###  plotting "newly quenched" quiescent fractions
        sp2.errorbar(lmasses_TOT_tstart[i_tau], f_newly_quenched_tend[i_tau],
                     color=colors[i_tau], lw=4, marker='o', ls=linestyles[i_tau%2], 
                     ms=8, ecolor='k', elinewidth=1.5, mew=2.,
                     label=taus_str[i_tau])


        ###  plotting matched quiescent fractions at z=1.05
        sp3.errorbar(lmasses_TOT_tstart[i_tau], fqu_tstart[i_tau],
                     color=colors[i_tau], lw=4, marker='o', ls=linestyles[i_tau%2], 
                     ms=8, ecolor='k', elinewidth=1.5, mew=2.,
                     label=taus_str[i_tau])


        ###  plotting predicted quiescent fractions at z=0.8
        sp4.errorbar(lmasses_TOT_tstart[i_tau], fqu_tend[i_tau],
                     color=colors[i_tau], lw=4, marker='o', ls=linestyles[i_tau%2], 
                     ms=8, ecolor='k', elinewidth=1.5, mew=2.,
                     label=taus_str[i_tau])

    sp1.legend(loc=1, fontsize=18, title=r'$\tau$ / Gyr')











'''





t0 = cosmo.age(0.8).value

lmass = numpy.array([numpy.log10(g.get_stellar_mass(t0)) for g in my_galaxies[1]])
sfqu = numpy.array([g.get_sfqu_flag(t0) for g in my_galaxies[1]])

digi_lmass = numpy.digitize(lmass, lmassbins)

ntot = []
nqu  = []

for i_lmass in range(len(lmassbins)):

    inds_tot = numpy.where(digi_lmass == i_lmass+1)[0]
    inds_qu  = numpy.where((digi_lmass == i_lmass+1) & (sfqu == 0))[0]

    ntot.append(len(inds_tot))
    nqu.append(len(inds_qu))

ntot = numpy.array(ntot)
nqu = numpy.array(nqu)






'''

















