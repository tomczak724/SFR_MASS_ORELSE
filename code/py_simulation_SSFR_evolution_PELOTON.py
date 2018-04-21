
import os
import time
import fsps
import glob
import numpy
import argparse
from os import path


###  setting up argument parser
if True:

    parser = argparse.ArgumentParser()

    parser.add_argument('-n1',
                        type=int,
                        help='Lower number for iteration IDs')

    parser.add_argument('-n2',
                        type=int,
                        help='Upper number for iteration IDs')

    parser.add_argument('-sfh',
                        help='Sets the Star-Formation History (SFH)',
                        choices=['single_exp', 'double_exp', 'single_exp_burst'])

    parser.add_argument('-tau_i',
                        type=float,
                        help='Value of tau for the initial portion of the SFH (Gyr)')

    parser.add_argument('-tau_f',
                        type=float,
                        help='Value of tau for the final portion of the SFH (Gyr)')

    parser.add_argument('-t_trans',
                        type=float,
                        help='Time after zstart=1.05 to transition from tau_i to tau_f (Gyr)')

    parser.add_argument('-t_burst',
                        type=float,
                        help='Time after zstart=1.05 to initiate the SF burst (Gyr)')

    parser.add_argument('-dt_burst',
                        type=float,
                        help='Duration of the SF burst (Gyr)')

    parser.add_argument('-scale_burst',
                        type=float,
                        help='Scale factor to multiply SFH during the SF burst')

    args = parser.parse_args()







###  establishing redshift-age relation for this cosmology
cosmo_data = numpy.loadtxt('/home/atomczak/PROJECTS/SFR-MASS_ORELSE/data/GRID_TIME_GRID_Z.dat')
#cosmo_data = numpy.loadtxt('../data/GRID_TIME_GRID_Z.dat')
universe_age = cosmo_data[:,0]
universe_redshifts = cosmo_data[:,1]




###  setting up FSPS model
print('\nGenerating FSPS stellar population model...')
t0 = time.time()

sps = fsps.StellarPopulation(zcontinuous=1)  # zcontinuous = 1 = interpolate SSPs to value given by "logzsol"
sps.params['logzsol'] = 0                 # Adopting solar metallicity
sps.params['sfh'] = 3                     # 1 = Tau model ; 3 = User input model
sps.params['imf_type'] = 1                # Chabrier (2003)
sps.params['add_igm_absorption'] = True   # Add IGM absorption
sps.params['dust_type'] = 2               # Calzetti+2000 attenuation curve

filternames = ['U', 'V', '2MASS_J']

tf = time.time()
print('done!  %.1f seconds\n' % (tf-t0))




def schechter_mf(xaxis, alpha, mstar, phistar, logphi=False):
    """
    #  DESCRIPTION:
    #    Returns the values for a Schechter mass function
    #    from a given mass-axis and Schechter parameters.
    #
    #  INPUTS:
    #    xaxis = input mass value(s)
    #    alpha = Schechter parameters
    #    mstar = Schechter parameters
    #    phistar = Schechter parameters
    #    logphi = Set True if phistar is Log(phistar)
    """
    if logphi:
        phistar = 10**phistar
    return numpy.log(10) * phistar * 10**((xaxis-mstar)*(1+alpha)) * numpy.exp(-10**(xaxis-mstar))


def dschechter_mf(xaxis, mstar, alpha1, phistar1, alpha2, phistar2, logphi=1):
    """
    #  DESCRIPTION:
    #    Returns the values for a Double Schechter mass function 
    #    from a given mass-axis and Schechter parameters.
    #
    #  INPUTS:
    #    xaxis  -------> input mass value(s)
    #    alpha[1,2] ---> Schechter parameters for [1,2] contribution
    #    mstar --------> Schechter parameter
    #    phistar[1,2] -> Schechter parameters for [1,2] contribution
    #    logphi -------> Set True if phistars are Log(phistar)
    """
    if logphi:
        phistar1 = 10**phistar1
        phistar2 = 10**phistar2
    sch1 = numpy.log(10) * phistar1 * 10**((xaxis-mstar)*(1+alpha1)) * numpy.exp(-10**(xaxis-mstar))
    sch2 = numpy.log(10) * phistar2 * 10**((xaxis-mstar)*(1+alpha2)) * numpy.exp(-10**(xaxis-mstar))
    return sch1 + sch2


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


def uvj_select( uv, vj, z ):
    """
    #  INPUT:
    #    uv  --->  Restframe (U-V) color (in AB system)  [array-like]
    #    vj  --->  Restframe (V-J) color (in AB system)  [array-like]
    #    z  ---->  Galaxy redshifts  [array-like]
    #
    #  OUTPUT: quiescent=0, star-forming=1
    #
    #  EXAMPLE:
    #            [ SF, SF, Qui, SF,  ....  Qui, Qui]
    #     array( [ 1,  1,  0,   1,   ....  0,   0 ] )
    #
    #  based on Whitaker+2011
    """

    floor = numpy.zeros(len(z))
    slope = numpy.zeros(len(z))
    intercept = numpy.zeros(len(z))

    floor[ numpy.where( (0.0<=z) & (z<1.5) ) ] = 1.3
    floor[ numpy.where( (1.5<=z) & (z<2.0) ) ] = 1.3
    floor[ numpy.where( (2.0<=z) ) ] = 1.2

    slope[ numpy.where( z<0.5 ) ] = 0.88
    slope[ numpy.where( 0.5<=z ) ] = 0.88

    intercept[ numpy.where( z<0.5 ) ] = 0.69
    intercept[ numpy.where( 0.5<=z ) ] = 0.59

    outer = numpy.zeros(len(z))
    outer[ numpy.where( (uv<slope*vj+intercept) | (uv<floor) )[0] ] = 1
    return outer.astype(int)




###  setting global variables
DT = 0.01
GRID_TIME = numpy.arange(DT, 12., DT)
GRID_Z = numpy.interp(GRID_TIME, 
                      universe_age, 
                      universe_redshifts)

OUTPUT_DIR = '/home/atomczak/PROJECTS/SFR-MASS_ORELSE/output'
#OUTPUT_DIR = '../data/simulations'





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
        return numpy.interp(tage, GRID_TIME, self.sfr)

    def get_ssfr(self, tage):
        '''
        Returns the SSFR at the requested age (in Gyr)
        '''
        return numpy.interp(tage, GRID_TIME, self.sfr/self.stellar_mass)

    def get_stellar_mass(self, tage):
        '''
        Returns the stellar mass at the requested age (in Gyr)
        '''
        return numpy.interp(tage, GRID_TIME, self.stellar_mass)

    def get_UV_color(self, tage):
        '''
        Returns the rest-frame (U-V) at the requested age (in Gyr)
        '''
        return numpy.interp(tage, GRID_TIME, self.UV_color)

    def get_VJ_color(self, tage):
        '''
        Returns the rest-frame (V-J) at the requested age (in Gyr)
        '''
        return numpy.interp(tage, GRID_TIME, self.VJ_color)

    def get_sfqu_flag(self, tage):
        '''
        Returns the SF/QU flag at the requested age (in Gyr)
        '''
        uv = self.get_UV_color(tage)
        vj = self.get_VJ_color(tage)
        z = numpy.interp(tage, universe_age, universe_redshifts)
        return uvj_select(uv, vj, numpy.array([z]))[0]










#######################################
###  Generating simulated galaxies  ###
###    1) Single exponential SFH    ###
#######################################

if args.sfh == 'single_exp':
    print('tau_i = %.4f' % args.tau_i)

    ###  dSchechter parameters for the SF-GSMF at 1.00<z<1.25 from Tomczak+2014
    lmstar = 10.44
    alpha1 = 0.53
    alpha2 = -1.44
    phistar1 = 10**(-2.98)
    phistar2 = 10**(-3.11)
    params = [lmstar, alpha1, phistar1, alpha2, phistar2]

    lmassax = numpy.linspace(8.5, 12., 1000)
    phiax = dschechter_mf(lmassax, *params)




    ###  selecting N random samples from the above GSMF
    N = args.n2 - args.n1

    phiax_cdf = numpy.array([phiax[:i].sum() for i in range(len(phiax))])
    phiax_cdf /= phiax_cdf.max()

    random_samples = numpy.random.rand(N)
    random_lmasses = numpy.interp(random_samples, phiax_cdf, lmassax)





    ###  applying SFRs based on Tomczak+2016 with +/- 0.2 dex scatter
    random_lsfrs = log_SFR_tomczak2016(1.05, random_lmasses, SF_galaxies=True)
    random_lsfrs += numpy.random.randn(N) * 0.2
    random_lssfrs = random_lsfrs - random_lmasses

    ###  Assigning exponential SFHs for random galaxies
    ###  Constraints are:
    ###    1) SFH prior to zstart=1.05 are exponential with fixed tau
    ###    2) The early part of this SFH must produce the instantaneous
    ###       SFR and M* already assigned based on the z=1.05 SFR-M* relation


    zstart = 1.05
    zend = 0.8

    tstart = numpy.interp(zstart, universe_redshifts[::-1], universe_age[::-1])
    tend = numpy.interp(zend, universe_redshifts[::-1], universe_age[::-1])




    i_gal = 0
    while i_gal < N:

        output_subdir = '%s/single_exp/%.4f' % (OUTPUT_DIR, args.tau_i)
        if not path.exists(output_subdir):
            os.mkdir(output_subdir)

        if path.exists('%s/galaxy%07i.txt' % (output_subdir, args.n1+i_gal)):
            i_gal += 1

        else:

            t0 = time.time()
            lmass_i = random_lmasses[i_gal]
            lsfr_i = random_lsfrs[i_gal]
            sfh0 = 10**lsfr_i * numpy.exp(-(GRID_TIME - tstart) / args.tau_i)



            ###  iterating though formation times until finding the specific tform
            ###  that reproduces the instantaneous SFR and M* for this galaxy
            lmass_array_tstarts = numpy.array([])
            for i_tform in range(len(GRID_TIME)):

                tform = GRID_TIME[i_tform]

                ###  for formation times after zstart=1.05 just add 0 to array
                if tform >= tstart:
                    lmass_array_tstarts = numpy.append(lmass_array_tstarts, 0.)
                else:
                    inds_truncate = numpy.where(GRID_TIME < tform)
                    sfh1 = sfh0 * 1.
                    sfh1[inds_truncate] = 0.
                    sps.set_tabular_sfh(GRID_TIME, sfh1*1.)

                    #lmass_tstart = numpy.interp(tstart, 10**(sps.log_age-9), numpy.log10(sps.stellar_mass))
                    sps.params['tage'] = tstart
                    lmass_tstart = numpy.log10(sps.stellar_mass)
                    lmass_array_tstarts = numpy.append(lmass_array_tstarts, lmass_tstart)


            ###  adopting best-fit formation time/redshift
            tform_best = numpy.interp(lmass_i, lmass_array_tstarts[::-1], GRID_TIME[::-1])
            inds_truncate = numpy.where(GRID_TIME < tform_best)
            sfh1 = sfh0 * 1.
            sfh1[inds_truncate] = 0.
            sps.set_tabular_sfh(GRID_TIME, sfh1*1.)


            ###  addind data to "galaxy" object
            g = galaxy()

            g.tau = args.tau_i * 1.
            g.dt = DT * 1.
            g.sfr = sfh1 * 1.

            ###  adding evolution of UVJ colors
            for i_time in range(len(GRID_TIME)):

                uvj_mags = sps.get_mags(tage=GRID_TIME[i_time], 
                                        redshift=0., 
                                        bands=filternames)

                uv = uvj_mags[0] - uvj_mags[1]
                vj = uvj_mags[1] - uvj_mags[2]
                sfqu = uvj_select(uv, vj, numpy.array([GRID_Z[i_time]]))[0]

                g.append_UV_color(uv)
                g.append_VJ_color(vj)
                g.append_sfqu_flag(sfqu)
                g.append_stellar_mass(sps.stellar_mass)

            tf = time.time()
            g.cpu_time = tf - t0



            ###  saving to ascii file
            outer = open('%s/galaxy%07i.txt' % (output_subdir, args.n1+i_gal), 'w')
            outer.write('# Data for sumlated galaxy %i\n' % (args.n1+i_gal))
            outer.write('# tau = %.4f\n' % args.tau_i)
            outer.write('# tform_best = %.4f\n' % tform_best)
            outer.write('# cpu_time = %.1f\n' % g.cpu_time)
            outer.write('# sfr stellar_mass UV_color VJ_color\n')
            for i in range(len(g.sfr)):
                outer.write(' %.4f' % g.sfr[i])
                outer.write(' %.4f' % g.stellar_mass[i])
                outer.write(' %.4f' % g.UV_color[i])
                outer.write(' %.4f' % g.VJ_color[i])
                outer.write('\n')
            outer.close()
            i_gal += 1













#######################################
###  Generating simulated galaxies  ###
###    2) Double exponential SFH    ###
#######################################

if args.sfh == 'double_exp':
    print('tau_i = %.4f' % args.tau_i)
    print('tau_f = %.4f' % args.tau_f)
    print('t_trans = %.4f' % args.t_trans)

    ###  dSchechter parameters for the SF-GSMF at 1.00<z<1.25 from Tomczak+2014
    lmstar = 10.44
    alpha1 = 0.53
    alpha2 = -1.44
    phistar1 = 10**(-2.98)
    phistar2 = 10**(-3.11)
    params = [lmstar, alpha1, phistar1, alpha2, phistar2]

    lmassax = numpy.linspace(8.5, 12., 1000)
    phiax = dschechter_mf(lmassax, *params)




    ###  selecting N random samples from the above GSMF
    N = args.n2 - args.n1

    phiax_cdf = numpy.array([phiax[:i].sum() for i in range(len(phiax))])
    phiax_cdf /= phiax_cdf.max()

    random_samples = numpy.random.rand(N)
    random_lmasses = numpy.interp(random_samples, phiax_cdf, lmassax)





    ###  applying SFRs based on Tomczak+2016 with +/- 0.2 dex scatter
    random_lsfrs = log_SFR_tomczak2016(1.05, random_lmasses, SF_galaxies=True)
    random_lsfrs += numpy.random.randn(N) * 0.2
    random_lssfrs = random_lsfrs - random_lmasses

    ###  Assigning exponential SFHs for random galaxies
    ###  Constraints are:
    ###    1) SFH prior to zstart=1.05 are exponential with fixed tau
    ###    2) The early part of this SFH must produce the instantaneous
    ###       SFR and M* already assigned based on the z=1.05 SFR-M* relation


    zstart = 1.05
    zend = 0.8

    tstart = numpy.interp(zstart, universe_redshifts[::-1], universe_age[::-1])
    tend = numpy.interp(zend, universe_redshifts[::-1], universe_age[::-1])




    i_gal = 0
    while i_gal < N:

        output_subdir = '%s/double_exp/%.4f_%.4f_%.4f' % (OUTPUT_DIR, args.tau_i, args.tau_f, args.t_trans)
        if not path.exists(output_subdir):
            os.mkdir(output_subdir)


        if path.exists('%s/galaxy%07i.txt' % (output_subdir, args.n1+i_gal)):
            i_gal += 1

        else:

            t0 = time.time()
            lmass_i = random_lmasses[i_gal]
            lsfr_i = random_lsfrs[i_gal]
            sfh0 = 10**lsfr_i * numpy.exp(-(GRID_TIME - tstart) / args.tau_i)



            ###  iterating though formation times until finding the specific tform
            ###  that reproduces the instantaneous SFR and M* for this galaxy
            lmass_array_tstarts = numpy.array([])
            for i_tform in range(len(GRID_TIME)):

                tform = GRID_TIME[i_tform]

                ###  for formation times after zstart=1.05 just add 0 to array
                if tform >= tstart:
                    lmass_array_tstarts = numpy.append(lmass_array_tstarts, 0.)
                else:
                    inds_truncate = numpy.where(GRID_TIME < tform)
                    sfh1 = sfh0 * 1.
                    sfh1[inds_truncate] = 0.
                    sps.set_tabular_sfh(GRID_TIME, sfh1*1.)

                    #lmass_tstart = numpy.interp(tstart, 10**(sps.log_age-9), numpy.log10(sps.stellar_mass))
                    sps.params['tage'] = tstart
                    lmass_tstart = numpy.log10(sps.stellar_mass)
                    lmass_array_tstarts = numpy.append(lmass_array_tstarts, lmass_tstart)


            ###  adopting best-fit formation time/redshift
            tform_best = numpy.interp(lmass_i, lmass_array_tstarts[::-1], GRID_TIME[::-1])
            inds_truncate = numpy.where(GRID_TIME < tform_best)
            sfh1 = sfh0 * 1.
            sfh1[inds_truncate] = 0.


            ###  applying secondary tau portion to SFH
            t_trans_cosmo = tstart+args.t_trans
            inds_tau_f = numpy.where(GRID_TIME > t_trans_cosmo)
            norm = sfh1[inds_tau_f][0]
            sfh1[inds_tau_f] = norm * numpy.exp(-(GRID_TIME[inds_tau_f]-t_trans_cosmo) / args.tau_f)


            sps.set_tabular_sfh(GRID_TIME, sfh1*1.)


            ###  addind data to "galaxy" object
            g = galaxy()

            g.tau = args.tau_i * 1.
            g.dt = DT * 1.
            g.sfr = sfh1 * 1.

            ###  adding evolution of UVJ colors
            for i_time in range(len(GRID_TIME)):

                uvj_mags = sps.get_mags(tage=GRID_TIME[i_time], 
                                        redshift=0., 
                                        bands=filternames)

                uv = uvj_mags[0] - uvj_mags[1]
                vj = uvj_mags[1] - uvj_mags[2]
                sfqu = uvj_select(uv, vj, numpy.array([GRID_Z[i_time]]))[0]

                g.append_UV_color(uv)
                g.append_VJ_color(vj)
                g.append_sfqu_flag(sfqu)
                g.append_stellar_mass(sps.stellar_mass)

            tf = time.time()
            g.cpu_time = tf - t0



            ###  saving to ascii file
            outer = open('%s/galaxy%07i.txt' % (output_subdir, args.n1+i_gal), 'w')
            outer.write('# Data for sumlated galaxy %i\n' % (args.n1+i_gal))
            outer.write('# tau_i = %.4f\n' % args.tau_i)
            outer.write('# tau_f = %.4f\n' % args.tau_f)
            outer.write('# t_trans = %.4f\n' % args.t_trans)
            outer.write('# cpu_time = %.1f\n' % g.cpu_time)
            outer.write('# sfr stellar_mass UV_color VJ_color\n')
            for i in range(len(g.sfr)):
                outer.write(' %.4f' % g.sfr[i])
                outer.write(' %.4f' % g.stellar_mass[i])
                outer.write(' %.4f' % g.UV_color[i])
                outer.write(' %.4f' % g.VJ_color[i])
                outer.write('\n')
            outer.close()
            i_gal += 1










###############################################
###  Generating simulated galaxies          ###
###    3) Single exponential with SF burst  ###
###############################################

if args.sfh == 'single_exp_burst':
    print('tau_i = %s' % args.tau_i)
    print('t_burst = %.4f' % args.t_burst)
    print('dt_burst = %.4f' % args.dt_burst)
    print('scale_burst = %.4f' % args.scale_burst)

    ###  dSchechter parameters for the SF-GSMF at 1.00<z<1.25 from Tomczak+2014
    lmstar = 10.44
    alpha1 = 0.53
    alpha2 = -1.44
    phistar1 = 10**(-2.98)
    phistar2 = 10**(-3.11)
    params = [lmstar, alpha1, phistar1, alpha2, phistar2]

    lmassax = numpy.linspace(8.5, 12., 1000)
    phiax = dschechter_mf(lmassax, *params)




    ###  selecting N random samples from the above GSMF
    N = args.n2 - args.n1

    phiax_cdf = numpy.array([phiax[:i].sum() for i in range(len(phiax))])
    phiax_cdf /= phiax_cdf.max()

    random_samples = numpy.random.rand(N)
    random_lmasses = numpy.interp(random_samples, phiax_cdf, lmassax)





    ###  applying SFRs based on Tomczak+2016 with +/- 0.2 dex scatter
    random_lsfrs = log_SFR_tomczak2016(1.05, random_lmasses, SF_galaxies=True)
    random_lsfrs += numpy.random.randn(N) * 0.2
    random_lssfrs = random_lsfrs - random_lmasses

    ###  Assigning exponential SFHs for random galaxies
    ###  Constraints are:
    ###    1) SFH prior to zstart=1.05 are exponential with fixed tau
    ###    2) The early part of this SFH must produce the instantaneous
    ###       SFR and M* already assigned based on the z=1.05 SFR-M* relation


    zstart = 1.05
    zend = 0.8

    tstart = numpy.interp(zstart, universe_redshifts[::-1], universe_age[::-1])
    tend = numpy.interp(zend, universe_redshifts[::-1], universe_age[::-1])




    i_gal = 0
    while i_gal < N:

        output_subdir = '%s/single_exp_burst/%.4f_%.4f_%.4f_%.4f' % (OUTPUT_DIR, args.tau_i, args.t_burst, args.dt_burst, args.scale_burst)
        if not path.exists(output_subdir):
            os.mkdir(output_subdir)


        if path.exists('%s/galaxy%07i.txt' % (output_subdir, args.n1+i_gal)):
            i_gal += 1

        else:

            t0 = time.time()
            lmass_i = random_lmasses[i_gal]
            lsfr_i = random_lsfrs[i_gal]
            sfh0 = 10**lsfr_i * numpy.exp(-(GRID_TIME - tstart) / args.tau_i)



            ###  iterating though formation times until finding the specific tform
            ###  that reproduces the instantaneous SFR and M* for this galaxy
            lmass_array_tstarts = numpy.array([])
            for i_tform in range(len(GRID_TIME)):

                tform = GRID_TIME[i_tform]

                ###  for formation times after zstart=1.05 just add 0 to array
                if tform >= tstart:
                    lmass_array_tstarts = numpy.append(lmass_array_tstarts, 0.)
                else:
                    inds_truncate = numpy.where(GRID_TIME < tform)
                    sfh1 = sfh0 * 1.
                    sfh1[inds_truncate] = 0.
                    sps.set_tabular_sfh(GRID_TIME, sfh1*1.)

                    #lmass_tstart = numpy.interp(tstart, 10**(sps.log_age-9), numpy.log10(sps.stellar_mass))
                    sps.params['tage'] = tstart
                    lmass_tstart = numpy.log10(sps.stellar_mass)
                    lmass_array_tstarts = numpy.append(lmass_array_tstarts, lmass_tstart)


            ###  adopting best-fit formation time/redshift
            tform_best = numpy.interp(lmass_i, lmass_array_tstarts[::-1], GRID_TIME[::-1])
            inds_truncate = numpy.where(GRID_TIME < tform_best)
            sfh1 = sfh0 * 1.
            sfh1[inds_truncate] = 0.


            ###  applying SF burst to SFH
            t_burst_start = tstart + args.t_burst
            t_burst_end = t_burst_start + args.dt_burst
            inds_burst = numpy.where((GRID_TIME > t_burst_start) & (GRID_TIME < t_burst_end))[0]
            sfh1[inds_burst] *= args.scale_burst


            sps.set_tabular_sfh(GRID_TIME, sfh1*1.)



            ###  addind data to "galaxy" object
            g = galaxy()

            g.tau = args.tau_i * 1.
            g.dt = DT * 1.
            g.sfr = sfh1 * 1.

            ###  adding evolution of UVJ colors
            for i_time in range(len(GRID_TIME)):

                uvj_mags = sps.get_mags(tage=GRID_TIME[i_time], 
                                        redshift=0., 
                                        bands=filternames)

                uv = uvj_mags[0] - uvj_mags[1]
                vj = uvj_mags[1] - uvj_mags[2]
                sfqu = uvj_select(uv, vj, numpy.array([GRID_Z[i_time]]))[0]

                g.append_UV_color(uv)
                g.append_VJ_color(vj)
                g.append_sfqu_flag(sfqu)
                g.append_stellar_mass(sps.stellar_mass)

            tf = time.time()
            g.cpu_time = tf - t0



            ###  saving to ascii file
            outer = open('%s/galaxy%07i.txt' % (output_subdir, args.n1+i_gal), 'w')
            outer.write('# Data for sumlated galaxy %i\n' % (args.n1+i_gal))
            outer.write('# tau = %.4f\n' % args.tau_i)
            outer.write('# t_burst = %.4f\n' % args.t_burst)
            outer.write('# dt_burst = %.4f\n' % args.dt_burst)
            outer.write('# scale_burst = %.4f\n' % args.scale_burst)
            outer.write('# cpu_time = %.1f\n' % g.cpu_time)
            outer.write('# sfr stellar_mass UV_color VJ_color\n')
            for i in range(len(g.sfr)):
                outer.write(' %.4f' % g.sfr[i])
                outer.write(' %.4f' % g.stellar_mass[i])
                outer.write(' %.4f' % g.UV_color[i])
                outer.write(' %.4f' % g.VJ_color[i])
                outer.write('\n')
            outer.close()
            i_gal += 1










































