
#import time
#import fsps
#import mypy
#import glob
#import numpy
import argparse
#from os import path


###  parsing input arguments
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
                    help='Duration of the SF burst (Gyr)')

parser.add_argument('-scale_burst',
                    help='Scale factor to multiply SFH during the SF burst')

args = parser.parse_args()




if args.sfh == 'single_exp':
    print 'tau_i = %s' % args.tau_i
if args.sfh == 'double_exp':
    print 'tau_i = %s' % args.tau_i
    print 'tau_f = %s' % args.tau_f
    print 't_trans = %s' % args.t_trans
if args.sfh == 'single_exp_burst':
    print 'tau_i = %s' % args.tau_i
    print 't_burst = %s' % args.t_burst
    print 'dt_burst = %s' % args.dt_burst
    print 'scale_burst = %s' % args.scale_burst

print 'n1 = %s' % args.n1
print 'n2 = %s' % args.n2
print 'n2-n1 = %s' % (args.n2 - args.n1)



'''
subparsers = parser.add_subparsers(title='sfh', 
                                   description='valid modes of usage', 
                                   dest='mode')

parser_gen.add_argument('--resolution', type=int, required=False, default=1,
            help="Resolution interpolation multiplier for NOAA topography")
'''


