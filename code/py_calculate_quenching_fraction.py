
import numpy



def get_e_ssfr(lssfr, e_lssfr):
    '''
    Returns the error on SSFR from the given log(SSFR) and e_log(SSFR)
    '''
    return e_lssfr * numpy.log(10) * 10**(lssfr)






###  All following quantites are for the low redshift bin (0.7 < z < 0.9)
###  lssfr_h ---> measured log(SSFR) of intermediate-M*, high-density galaxies
###  lssfr_li --> measured log(SSFR) of intermediate-M*, low- and intermediate-density galaxies
###  lssfr_Q ---> simulated log(SSFR) of tau=0.5 Gyr galaxies

lssfr_h  = -9.6098
lssfr_li = -9.3100
lssfr_Q  = -9.8550

e_lssfr_h  = 0.0934
e_lssfr_li = 0.0657


ssfr_h = 10**lssfr_h
ssfr_li = 10**lssfr_li
ssfr_Q = 10**lssfr_Q

e_ssfr_h  = get_e_ssfr(lssfr_h, e_lssfr_h)
e_ssfr_li = get_e_ssfr(lssfr_li, e_lssfr_li)



fquenching = (ssfr_h - ssfr_li) / (ssfr_Q - ssfr_li)

e_fquenching2 = (e_ssfr_h / (ssfr_Q - ssfr_li))**2 \
                + e_ssfr_li**2 * (ssfr_h - ssfr_Q)**2 / (ssfr_Q - ssfr_li)**4



print('\nfquenching = %.1f +/- %.1f %%\n' % (fquenching*100, e_fquenching2**0.5*100))





























