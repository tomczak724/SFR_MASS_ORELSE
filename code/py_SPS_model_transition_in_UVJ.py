
import os
import time
import fsps
import mypy
import numpy
import subprocess
from astropy import cosmology
from matplotlib import pyplot
from astropy.io import fits, ascii
from matplotlib.backends.backend_pdf import PdfPages

cosmo = cosmology.FlatLambdaCDM(H0=70., Om0=0.3)


###  establishing redshift-age relation for this cosmology
universe_redshifts = numpy.linspace(30, 0, 2000)
universe_age = cosmo.age(universe_redshifts)



class galaxy(object):
    def __init__(self):
        self.tau = 0.       # tau for exp-declining SFH
        self.dt = 0.
        self.t_trans = 0.
        self.time = numpy.array([])
        self.age_Gyr = numpy.array([])
        self.sfr = numpy.array([])
        self.stellar_mass = numpy.array([])
        self.UV_color = numpy.array([])
        self.VJ_color = numpy.array([])
        self.sfqu_flag = numpy.array([])

    def append_time(self, t):
        self.time = numpy.append(self.time, t)

    def append_sfr(self, s):
        self.sfr = numpy.append(self.sfr, s)

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
        return numpy.interp(tage, self.age_Gyr, self.sfr)

    def get_ssfr(self, tage):
        '''
        Returns the SSFR at the requested age (in Gyr)
        '''
        return numpy.interp(tage, self.age_Gyr, self.sfr/self.stellar_mass)

    def get_stellar_mass(self, tage):
        '''
        Returns the stellar mass at the requested age (in Gyr)
        '''
        return numpy.interp(tage, self.age_Gyr, self.stellar_mass)

    def get_UV_color(self, tage):
        '''
        Returns the rest-frame (U-V) at the requested age (in Gyr)
        '''
        return numpy.interp(tage, self.age_Gyr, self.UV_color)

    def get_VJ_color(self, tage):
        '''
        Returns the rest-frame (V-J) at the requested age (in Gyr)
        '''
        return numpy.interp(tage, self.age_Gyr, self.VJ_color)

    def get_sfqu_flag(self, tage):
        '''
        Returns the SF/QU flag at the requested age (in Gyr)
        '''
        uv = self.get_UV_color(tage)
        vj = self.get_VJ_color(tage)
        z = numpy.interp(tage, universe_age, universe_redshifts)
        return mypy.uvj_select(uv, vj, [z])[0]

    def smoothed_sfr(self, dt=0.1):
        '''
        Returns the SFR array smoothed by a boxcar kernel
        of the given width dt (in Gyr)
        '''
        k = convolution.Box1DKernel(width=dt)
        return convolution.convolve(self.sfr, k)


class galaxy2(object):
    def __init__(self):
        self.tau1 = 0.       # tau for initail exp-declining SFH
        self.tau2 = 0.       # tau for exp-decline after SF "burst"
        self.tburst = 0.     # time at which SF "burst" occurs (Gyr)
        self.dt = 0.         # discretized time unit (Gyr)

        self.time = numpy.array([])
        self.sfr = numpy.array([])
        self.stellar_mass = numpy.array([])
        self.UV_color = numpy.array([])
        self.VJ_color = numpy.array([])
        self.sfqu_flag = numpy.array([])

    def append_time(self, t):
        self.time = numpy.append(self.time, t)

    def append_sfr(self, s):
        self.sfr = numpy.append(self.sfr, s)

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
        Returns the SFR at the requested time (in Gyr)
        '''
        return numpy.interp(tage, self.time, self.sfr)

    def get_ssfr(self, tage):
        '''
        Returns the SSFR at the requested time (in Gyr)
        '''
        return numpy.interp(tage, self.time, self.sfr/self.stellar_mass)

    def get_stellar_mass(self, tage):
        '''
        Returns the stellar mass at the requested time (in Gyr)
        '''
        return numpy.interp(tage, self.time, self.stellar_mass)

    def get_UV_color(self, tage):
        '''
        Returns the rest-frame (U-V) at the requested time (in Gyr)
        '''
        return numpy.interp(tage, self.time, self.UV_color)

    def get_VJ_color(self, tage):
        '''
        Returns the rest-frame (V-J) at the requested time (in Gyr)
        '''
        return numpy.interp(tage, self.time, self.VJ_color)

    def get_sfqu_flag(self, tage):
        '''
        Returns the SF/QU flag at the requested time (in Gyr)
        '''
        uv = self.get_UV_color(tage)
        vj = self.get_VJ_color(tage)
        z = numpy.interp(tage, universe_age, universe_redshifts)
        return mypy.uvj_select(uv, vj, [z])[0]

    def smoothed_sfr(self, tsmooth=0.1):
        '''
        Returns the SFR array smoothed by a boxcar kernel
        of the given width dt (in Gyr)
        '''
        smooth_length = int(tsmooth / dt + 0.5)
        k = convolution.Box1DKernel(width=smooth_length)
        return convolution.convolve(self.sfr, k)






print '\nGenerating FSPS stellar population models...'
t0 = time.time()

sps = fsps.StellarPopulation(zcontinuous=1)

sps.params['sfh'] = 3                     # 1 = Tau model ; 3 = User input model
sps.params['imf_type'] = 1                # Chabrier (2003)
sps.params['add_igm_absorption'] = True   # Add IGM absorption
sps.params['dust_type'] = 2               # Calzetti+2000 attenuation curve

filternames = ['U', 'V', '2MASS_J']

tf = time.time()
print 'done!  %.1f seconds\n' % (tf-t0)





###  Plotting cartoon figure to illustrate the const+exp SFH
if False:

    fig = pyplot.figure()

    sp = fig.add_subplot(111)
    sp.set_xlabel('time', size=25)
    sp.set_ylabel('log( SFR )', size=25)

    sp.set_xticks([])
    sp.set_yticks([])

    sp.axis([3, 8, -1.2, 0.7])

    sfh = numpy.ones(N_AGE)
    sfh = numpy.zeros(N_AGE)

    dt = 2.
    tau = 1.

    t_trans = cosmo.age(0.9).value - dt
    inds_trans = numpy.where(GRID_AGE > t_trans)[0]
    sfh[inds_trans] =  numpy.exp(-(GRID_AGE[inds_trans]-t_trans) / tau)


    ###  plotting SFH
    sp.plot(GRID_AGE, numpy.log10(sfh), color='b', lw=4)

    ###  plotting vline for z=0.9
    sp.axvline(cosmo.age(0.9).value, color='gray', ls=':')
    t1 = sp.text(cosmo.age(0.9).value * 1.02,
                 numpy.log10(sfh).max() - 0.2,
                 'z = 0.9',
                 horizontalalignment='left',
                 verticalalignment='top')


    x0 = cosmo.age(0.9).value - dt/2.
    y0 = numpy.log10(sfh).max() + 0.2
    dx = dt / 2.

    sp.errorbar(x0, y0, xerr=dx, ecolor='k', elinewidth=3, capsize=7)
    t2 = sp.text(x0, 
                 y0, 
                 u'$\Delta$t',
                 horizontalalignment='center',
                 verticalalignment='bottom')


    y1 = numpy.interp(x0, GRID_AGE, numpy.log10(sfh))
    t3 = sp.text(x0,
                 y1,
                 r'$\propto$ exp( $-$t / $\tau$ )   ',
                 horizontalalignment='right',
                 verticalalignment='top')

    fig.savefig('../figures/SFH_constant+exponential.png')
    pyplot.close()








###  setting global variables
DT = 0.01
GRID_TIME = numpy.arange(DT, 12., DT)
GRID_Z = numpy.interp(GRID_TIME, 
                      universe_age, 
                      universe_redshifts)






tau0 = 2.
norm0 = 10.
sfh0 = norm0 * numpy.exp(-GRID_TIME / tau0)
sps.set_tabular_sfh(GRID_TIME, sfh0*1.)


my_galaxies = []
g = galaxy2()

for i_time in range(len(GRID_TIME)):

    uvj_mags = sps.get_mags(tage=GRID_TIME[i_time], 
                            redshift=0., 
                            bands=filternames)

    uv = uvj_mags[0] - uvj_mags[1]
    vj = uvj_mags[1] - uvj_mags[2]
    sfqu = mypy.uvj_select(uv, vj, [GRID_Z[i_time]])[0]

    ###  ONCE THE GALAXY HAS QUIESCENT UVJ COLORS STOP!!!
    if sfqu == 0:
        break


    g.append_time(GRID_TIME[i_time])
    g.append_sfr(sfh0[i_time])
    g.append_UV_color(uv)
    g.append_VJ_color(vj)
    g.append_sfqu_flag(sfqu)
    g.append_stellar_mass(sps.stellar_mass)





sim_lmass = numpy.random.rand(1000) * 0.7 + 9.4













my_galaxies = []

###  setting up progressbar
n_iterations = len(GRID_DT) * len(GRID_TAU)
n_iterations = 29 * len(GRID_TAU)
widgets = ['  Running: ', progressbar.Percentage(), 
           ' ', progressbar.Bar(marker=progressbar.RotatingMarker()), 
           ' ', progressbar.ETA(), 
           ' ', progressbar.FileTransferSpeed()]
pbar = progressbar.ProgressBar(widgets=widgets, maxval=n_iterations).start()

counter = -1
t0 = time.time()

for i_tau in range(len(GRID_TAU)):

    normalization = 1.0

    sfh = numpy.ones(N_AGE) * normalization
    sfh = numpy.zeros(N_AGE) * normalization


    tau = GRID_TAU[i_tau]
    GRID_DT = numpy.arange(tau/10., 3*tau, tau/10.)

    for i_dt in range(len(GRID_DT)):

        counter += 1
        pbar.update(counter)



        ###  identifying time-range for exponential phase
        t_trans = cosmo.age(0.9).value - GRID_DT[i_dt]
        inds_trans = numpy.where(GRID_AGE > t_trans)[0]

        sfh[inds_trans] = normalization * numpy.exp(-(GRID_AGE[inds_trans]-t_trans) / GRID_TAU[i_tau])
        sps.set_tabular_sfh(GRID_AGE, sfh*1.)

        my_galaxies.append(galaxy())

        my_galaxies[-1].tau = GRID_TAU[i_tau]
        my_galaxies[-1].dt = GRID_DT[i_dt]
        my_galaxies[-1].t_trans = t_trans
        my_galaxies[-1].age_Gyr = GRID_AGE
        my_galaxies[-1].redshifts = GRID_Z
        my_galaxies[-1].sfr = sfh * 1.




        for i_age in range(len(GRID_AGE)):

            uvj_mags = sps.get_mags(tage=GRID_AGE[i_age], 
                                    redshift=0., 
                                    bands=filternames)

            uv = uvj_mags[0] - uvj_mags[1]
            vj = uvj_mags[1] - uvj_mags[2]

            my_galaxies[-1].append_UV_color(uv)
            my_galaxies[-1].append_VJ_color(vj)

            sfqu = mypy.uvj_select(uv, vj, [GRID_Z[i_age]])[0]
            my_galaxies[-1].append_sfqu_flag(sfqu)

            my_galaxies[-1].append_stellar_mass(sps.stellar_mass)

tf = time.time()
print '\n\n%.1f minutes' % ((tf-t0)/60.)
















my_galaxies = []

###  time grid for the age of the universe
GRID_TIME = numpy.arange(0.001, 12., 0.001)




###  setting up progressbar
n_iterations = len(GRID_TAU)
widgets = ['  Running: ', progressbar.Percentage(), 
           ' ', progressbar.Bar(marker=progressbar.RotatingMarker()), 
           ' ', progressbar.ETA(), 
           ' ', progressbar.FileTransferSpeed()]
pbar = progressbar.ProgressBar(widgets=widgets, maxval=n_iterations).start()



for i_tau in range(len(GRID_TAU)):

    pbar.update(i_tau)

    tau = GRID_TAU[i_tau]
    dt = tau / 10.

    sfh = numpy.exp(-GRID_TIME / tau)
    sps.set_tabular_sfh(GRID_TIME, sfh*1.)


    my_galaxies.append(galaxy())

    my_galaxies[-1].tau = tau
    my_galaxies[-1].dt = dt


    for i_t in range(len(GRID_TIME)):

        uvj_mags = sps.get_mags(tage=GRID_TIME[i_t], 
                                redshift=0., 
                                bands=filternames)

        uv = uvj_mags[0] - uvj_mags[1]
        vj = uvj_mags[1] - uvj_mags[2]
        sfqu = mypy.uvj_select(uv, vj, [0.9])[0]
        if sfqu == 0:
            break

        my_galaxies[-1].append_time(GRID_TIME[i_t])
        my_galaxies[-1].append_sfr(sfh[i_t])
        my_galaxies[-1].append_UV_color(uv)
        my_galaxies[-1].append_VJ_color(vj)
        my_galaxies[-1].append_sfqu_flag(sfqu)
        my_galaxies[-1].append_stellar_mass(sps.stellar_mass)









###  plotting SSFR vs. tform in bins of tau

fig = pyplot.figure()

sp = fig.add_subplot(111)

sp.minorticks_on()
sp.set_xlabel('$\Delta$t$_{form}$')
sp.set_ylabel('log( SSFR$_{\,z=0.9}$ / [yr$^{-1}$] )')

fig.subplots_adjust(top=0.89)

top_ax = sp.twiny()
top_ax.set_xlabel('formation redshift')



dt_form_axis = numpy.arange(0.1, 5.1, 0.1)

for i_gal in range(len(my_galaxies)):

    g = my_galaxies[i_gal]

    ###  selecting time-range over which galaxy is UVJ-SF
    inds_sf = numpy.where(dt_form_axis <= g.time.max())[0]


    ssfr0 = g.sfr / g.stellar_mass
    ssfr1 = numpy.interp(dt_form_axis, g.time, ssfr0)



    sp.plot(dt_form_axis[inds_sf], 
            numpy.log10(ssfr1[inds_sf]),
            marker='o')


























