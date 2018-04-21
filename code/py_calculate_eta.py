
import mypy
import numpy
import progressbar
from astropy.io import ascii
from astropy import units, constants, cosmology, coordinates

cosmo = cosmology.FlatLambdaCDM(H0=70., Om0=0.3)


class combined_catalog(object):
    def __init__(self, header, n_galaxies):

        self.time_stamp = time.ctime()
        self.header = header

        self.n_galaxies = n_galaxies



galaxies_catalog = pickle.load(open('../data/table_all_galaxies.pickle', 'rb'))
substructure_catalog = ascii.read('../data/ORELSE_substructures.dat')

substructure_centroids = coordinates.SkyCoord(substructure_catalog['RA'], substructure_catalog['Dec'], unit=units.deg, frame='fk5')
R200 = 3**0.5 * substructure_catalog['sigma_1Mpc']*(units.km/units.s) / \
       (10 * cosmo.H(substructure_catalog['z']))
substructure_velocities = mypy.vel_recess(substructure_catalog['z'])






eta = numpy.zeros(galaxies_catalog.n_galaxies) * numpy.nan




###  Initializing progress bar
widgets = ['  Running: ', progressbar.Percentage(), 
           ' ', progressbar.Bar(marker=progressbar.RotatingMarker()), 
           ' ', progressbar.ETA(), 
           ' ', progressbar.FileTransferSpeed()]
pbar = progressbar.ProgressBar(widgets=widgets, maxval=galaxies_catalog.n_galaxies+1).start()


for i_gal in range(galaxies_catalog.n_galaxies):

    pbar.update(i_gal)

    ###  calculate Rproj to every ORELSE substructure (that's right, every single one)
    galaxy_centroid = coordinates.SkyCoord(galaxies_catalog.ra[i_gal], galaxies_catalog.dec[i_gal], unit=units.degree, frame='fk5')

    separations = galaxy_centroid.separation(substructure_centroids)
    separations = separations.to(units.arcmin)

    Rproj = separations * cosmo.kpc_proper_per_arcmin(substructure_catalog['z'])
    Rproj = Rproj.to(units.Mpc)


    ###  calculating recessional velocities
    galaxy_velocity = mypy.vel_recess(galaxies_catalog.z[i_gal])
    delta_velocities = abs(galaxy_velocity - substructure_velocities)


    ###  calculating all possible etas
    eta_all = (Rproj / R200).value * (delta_velocities / substructure_catalog['sigma_1Mpc'])



    ###  1) Find all substructures within Rproj < 5 Mpc and adopt smallest eta
    inds_5Mpc = numpy.where(Rproj < 5*units.Mpc)[0]

    if len(inds_5Mpc) > 0:
        eta[i_gal] = eta_all[inds_5Mpc].min()









