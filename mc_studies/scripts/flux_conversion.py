import math
from scipy.integrate import quad
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.font_manager as fm
import matplotlib.ticker as mticker
prop_font = fm.FontProperties(fname='/home/hunt-stokes/LiberationSerif-Regular.ttf',size=28)

matplotlib.rcParams['xtick.major.size'] = 10
matplotlib.rcParams['xtick.major.width'] = 2
matplotlib.rcParams['xtick.minor.size'] = 5
matplotlib.rcParams['xtick.minor.width'] = 1
matplotlib.rcParams['ytick.major.size'] = 10
matplotlib.rcParams['ytick.major.width'] = 2
matplotlib.rcParams['ytick.minor.size'] = 5
matplotlib.rcParams['ytick.minor.width'] = 1
matplotlib.rcParams['axes.linewidth'] = 2 #set the value globally
matplotlib.rcParams['figure.facecolor'] = 'white'
matplotlib.rcParams['figure.figsize'] = 18, 10
matplotlib.rcParams['xtick.major.pad']='12'
matplotlib.rcParams['ytick.major.pad']='12'
#matplotlib.rcParams.update({'mathtext.default':  prop_font })

############################################################################
############################################################################
####
####  Section to calculate the cross sections with radiative corrections
####
############################################################################
############################################################################

def integrand(x):
  return float(math.log(1-x)/x)
def integrand2(x):
  return float(math.log(x)/(1-x))

def xSection(T):

  rhoNC = 1.0126
  rhoNC_err = 0.0016

  # Term multiplying the parenthesis
  term1 = 2*GF2*me / Pi

  x = math.sqrt( 1. + ( 2.*me/T ))
  x2 = x*x
  I = (1./6.) * ( (1./3.) + (3.-x2)*( (0.5*x*math.log((x+1.)/(x-1.))) -1. ))

  k = 0.
  k_err = 0.
  gL = 0.
  gR = 0.

  if type == 'nue':
    k = 0.9791 + (0.0097*I)
    k_err = 0.0025

    gL = (rhoNC * ( 0.5 - k*sin2_thetaW )) - 1
    gR = -rhoNC * k * sin2_thetaW

  if type == 'numu':
    k = 0.9970 + (0.00037*I)
    k_err = 0.0025

    gL = rhoNC * ( 0.5 - k*sin2_thetaW )
    gR = -rhoNC * k * sin2_thetaW

  gL2 = gL * gL
  gR2 = gR * gR

  z = T/Enu
  E = T + me

  l = math.sqrt((E*E) - (me*me))
  beta = l/E

  zmax = Tmax/Enu
  m1zmax = 1. - zmax
  Lz = quad(integrand, 0, zmax)
  L1mz = quad(integrand2, 0, m1zmax)

  Emax = Tmax + me
  betamax = math.sqrt((Emax*Emax) - (me*me))/Emax
  Lbeta = quad(integrand, 0, betamax)

  fminus = ( ((E/l)*math.log((E+l)/me)) -1. ) * ( 2.*math.log(1.-z-(me/(E+l))) - math.log(1.-z) - 0.5*math.log(z) - (5./12.) ) + 0.5*Lz[0] - 0.5*Lbeta[0] - 0.5*math.log(1.-z)*math.log(1.-z) - ((11./12.)+(z/2.))*math.log(1.-z) + z*(math.log(z) + 0.5*math.log(2.*Enu/me)) - ((31./18.)+(1./12.)*math.log(z))*beta - (11./12.)*z + (z*z/24.)

  oneminusz2fplus = ( (E/l)*math.log((E+l)/me) -1 ) * ( (1-z)*(1-z)*( 2*math.log(1-z-(me/(E+l))) - math.log(1-z) - (math.log(z)/2) - (2/3) ) - 0.5*(z*z*math.log(z) + 1 - z))  - 0.5*(1-z)*(1-z)*( math.log(1-z)*math.log(1-z) + beta*(L1mz[0] - math.log(z)*math.log(1-z))) + math.log(1-z) * ( 0.5*z*z*math.log(z) + (1-z)*(2*z-0.5)/3) - 0.5*z*z*L1mz[0] - (z*(1-2*z)*math.log(z)/3) - (z*(1-z)/6) - (beta/12)*(math.log(z) + ((1-z)*(115-109*z)/6))

  fpm = ((((E/l)*math.log((E+l)/me)))-1)*2*math.log( 1 -z - (me/(E+l)) )

  term2 = gL2*( 1 + ((alpha/Pi)*fminus))
  term3 = (gR2*(1-z)*(1-z)) + (gR2*(alpha/Pi)*oneminusz2fplus)
  term4 = -gR*gL*me*(z/Enu)*(1+((alpha/Pi)*fpm))

  dsigmadT = term1 * ( term2 + term3 + term4 ) *hbar*hbar*cl*cl

  return dsigmadT

############################################################################
############################################################################
####
####  Main Code
####
############################################################################
############################################################################

######################## Definition of Variables ###########################
# Physics constants
Pi = math.pi              # Pi
nA = 6.02214076E+23       # mol^-1 : Avogadro's Number
hbar = 6.582119569e-22    # MeV s : h-bar
alpha = 7.29735e-3           # Fine Structure Constant
GF = 1.1663787e-11        # MeV^-2 : Fermi Coupling Constant
GF2 = GF*GF               # MeV^-4 : Fermi Coupling Constant Squared
sin2_thetaW = 0.23116     # Weak Mixing Angle

cl = 29979245800 ## cm/s

me = 0.511  # MeV c^-2 : Electron Mass

# LAB+PPO characteristics
rho_LAB = 0.860    # g/cm^-3 : LAB density
m_weight = 235.0   # g/mol : LAB molecular weight
av_radius = 600.53 # cm : AV radius

# Oscillation Parameters from:
# PDG 2019 - M. Tanabashi et al. (Particle Data Group), Phys. Rev. D 98, 030001 (2018) and 2019 update.
sin2_12 = 0.307
sin2_12_err = 0.013

sin2_13 = 2.18e-2
sin2_13_err = 0.07e-2

sin2_23_i = 0.536      # Inverted order
sin2_23_i_errp = 0.023 # Inverted order, positive error bar
sin2_23_i_errn = 0.028 # Inverted order, negative error bar

sin2_23_n1 = 0.512      # Normal order octant I
sin2_23_n1_errp = 0.019 # Normal order octant I, positive error bar
sin2_23_n1_errn = 0.022 # Normal order octant I, negative error bar

sin2_23_n2 = 0.542      # Normal order octant II
sin2_23_n2_errp = 0.019 # Normal order octant II, positive error bar
sin2_23_n2_errn = 0.022 # Normal order octant II, negative error bar

dm2_12 = 7.53e-5     # eV^2
dm2_12_err = 0.18e-5 # eV^2

dm2_32_i = -2.55e-3     # eV^2 (Inverted order)
dm2_32_i_err = 0.04e-3  # eV^2 (Inverted order)
dm2_32_n = 2.444e-3     # eV^2 (Normal order)
dm2_32_n_err = 0.034e-3 # eV^2 (Normal order)


############################ Solar Neutrinos ###############################

E_7Be = [15]  # MeV : Energy of the 7Be Solar Neutrinos
BR_7Be = [1] # Branching ratios [0.897, 0.103]
flux_7Be = 5.46e6       # cm^-2 s^-1 : Flux at Earth of the 7Be Solar Neutrinos

############ Number of electron targets in the Scintillator ################

# Atoms in the LAB and PPO molecules - C, H, N and O
# Number of electrons in each:
atC = 6
atH = 1
atN = 7
atO = 8

# Molecules, fraction by mass, and cocktail content (99.9% LAB, 0.056% PPO for 0.5g/L in the PF)
molecules = [16*atC+26*atH,17*atC+28*atH,18*atC+30*atH,19*atC+32*atH,15*atC+24*atH,15*atC+11*atH+atN+atO]
content = [0.204, 0.432, 0.334, 0.018, 0.012, 1]
cocktail = [ 0.999425, 0.999425, 0.999425, 0.999425, 0.999425, 0.000575]

n_cocktail = 0

for i in range(len(molecules)):
  n_cocktail += cocktail[i]*content[i]*molecules[i]

print(' = ', n_cocktail)

n_cocktail = 0.999425 * 131 + 0.000575 * 116
print(' = ', n_cocktail)
mass = 365e6 # g

# Finally calculate the number of electron targets inside the AV volume
n_target = n_cocktail * nA * mass * (1./m_weight)

# Now want the number per 100t of LAB
#n_target = n_target * 1.E+8 / ( rho_LAB * vol )

print('Number of electron targets = ', n_target)

#exit()

######################### Energy Spectra from Bahcall ########################
energyB = [0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1,1.02,1.04,1.06,1.08,1.1,1.12,1.14,1.16,1.18,1.2,1.22,1.24,1.26,1.28,1.3,1.32,1.34,1.36,1.38,1.4,1.42,1.44,1.46,1.48,1.5,1.52,1.54,1.56,1.58,1.6,1.62,1.64,1.66,1.68,1.7,1.72,1.74,1.76,1.78,1.8,1.82,1.84,1.86,1.88,1.9,1.92,1.94,1.96,1.98,2,2.02,2.04,2.06,2.08,2.1,2.12,2.14,2.16,2.18,2.2,2.22,2.24,2.26,2.28,2.3,2.32,2.34,2.36,2.38,2.4,2.42,2.44,2.46,2.48,2.5,2.52,2.54,2.56,2.58,2.6,2.62,2.64,2.66,2.68,2.7,2.72,2.74,2.76,2.78,2.8,2.82,2.84,2.86,2.88,2.9,2.92,2.94,2.96,2.98,3,3.02,3.04,3.06,3.08,3.1,3.12,3.14,3.16,3.18,3.2,3.22,3.24,3.26,3.28,3.3,3.32,3.34,3.36,3.38,3.4,3.42,3.44,3.46,3.48,3.5,3.52,3.54,3.56,3.58,3.6,3.62,3.64,3.66,3.68,3.7,3.72,3.74,3.76,3.78,3.8,3.82,3.84,3.86,3.88,3.9,3.92,3.94,3.96,3.98,4,4.02,4.04,4.06,4.08,4.1,4.12,4.14,4.16,4.18,4.2,4.22,4.24,4.26,4.28,4.3,4.32,4.34,4.36,4.38,4.4,4.42,4.44,4.46,4.48,4.5,4.52,4.54,4.56,4.58,4.6,4.62,4.64,4.66,4.68,4.7,4.72,4.74,4.76,4.78,4.8,4.82,4.84,4.86,4.88,4.9,4.92,4.94,4.96,4.98,5,5.02,5.04,5.06,5.08,5.1,5.12,5.14,5.16,5.18,5.2,5.22,5.24,5.26,5.28,5.3,5.32,5.34,5.36,5.38,5.4,5.42,5.44,5.46,5.48,5.5,5.52,5.54,5.56,5.58,5.6,5.62,5.64,5.66,5.68,5.7,5.72,5.74,5.76,5.78,5.8,5.82,5.84,5.86,5.88,5.9,5.92,5.94,5.96,5.98,6,6.02,6.04,6.06,6.08,6.1,6.12,6.14,6.16,6.18,6.2,6.22,6.24,6.26,6.28,6.3,6.32,6.34,6.36,6.38,6.4,6.42,6.44,6.46,6.48,6.5,6.52,6.54,6.56,6.58,6.6,6.62,6.64,6.66,6.68,6.7,6.72,6.74,6.76,6.78,6.8,6.82,6.84,6.86,6.88,6.9,6.92,6.94,6.96,6.98,7,7.02,7.04,7.06,7.08,7.1,7.12,7.14,7.16,7.18,7.2,7.22,7.24,7.26,7.28,7.3,7.32,7.34,7.36,7.38,7.4,7.42,7.44,7.46,7.48,7.5,7.52,7.54,7.56,7.58,7.6,7.62,7.64,7.66,7.68,7.7,7.72,7.74,7.76,7.78,7.8,7.82,7.84,7.86,7.88,7.9,7.92,7.94,7.96,7.98,8,8.02,8.04,8.06,8.08,8.1,8.12,8.14,8.16,8.18,8.2,8.22,8.24,8.26,8.28,8.3,8.32,8.34,8.36,8.38,8.4,8.42,8.44,8.46,8.48,8.5,8.52,8.54,8.56,8.58,8.6,8.62,8.64,8.66,8.68,8.7,8.72,8.74,8.76,8.78,8.8,8.82,8.84,8.86,8.88,8.9,8.92,8.94,8.96,8.98,9,9.02,9.04,9.06,9.08,9.1,9.12,9.14,9.16,9.18,9.2,9.22,9.24,9.26,9.28,9.3,9.32,9.34,9.36,9.38,9.4,9.42,9.44,9.46,9.48,9.5,9.52,9.54,9.56,9.58,9.6,9.62,9.64,9.66,9.68,9.7,9.72,9.74,9.76,9.78,9.8,9.82,9.84,9.86,9.88,9.9,9.92,9.94,9.96,9.98,10,10.02,10.04,10.06,10.08,10.1,10.12,10.14,10.16,10.18,10.2,10.22,10.24,10.26,10.28,10.3,10.32,10.34,10.36,10.38,10.4,10.42,10.44,10.46,10.48,10.5,10.52,10.54,10.56,10.58,10.6,10.62,10.64,10.66,10.68,10.7,10.72,10.74,10.76,10.78,10.8,10.82,10.84,10.86,10.88,10.9,10.92,10.94,10.96,10.98,11,11.02,11.04,11.06,11.08,11.1,11.12,11.14,11.16,11.18,11.2,11.22,11.24,11.26,11.28,11.3,11.32,11.34,11.36,11.38,11.4,11.42,11.44,11.46,11.48,11.5,11.52,11.54,11.56,11.58,11.6,11.62,11.64,11.66,11.68,11.7,11.72,11.74,11.76,11.78,11.8,11.82,11.84,11.86,11.88,11.9,11.92,11.94,11.96,11.98,12,12.02,12.04,12.06,12.08,12.1,12.12,12.14,12.16,12.18,12.2,12.22,12.24,12.26,12.28,12.3,12.32,12.34,12.36,12.38,12.4,12.42,12.44,12.46,12.48,12.5,12.52,12.54,12.56,12.58,12.6,12.62,12.64,12.66,12.68,12.7,12.72,12.74,12.76,12.78,12.8,12.82,12.84,12.86,12.88,12.9,12.92,12.94,12.96,12.98,13,13.02,13.04,13.06,13.08,13.1,13.12,13.14,13.16,13.18,13.2,13.22,13.24,13.26,13.28,13.3,13.32,13.34,13.36,13.38,13.4,13.42,13.44,13.46,13.48,13.5,13.52,13.54,13.56,13.58,13.6,13.62,13.64,13.66,13.68,13.7,13.72,13.74,13.76,13.78,13.8,13.82,13.84,13.86,13.88,13.9,13.92,13.94,13.96,13.98,14,14.02,14.04,14.06,14.08,14.1,14.12,14.14,14.16,14.18,14.2,14.22,14.24,14.26,14.28,14.3,14.32,14.34,14.36,14.38,14.4,14.42,14.44,14.46,14.48,14.5,14.52,14.54,14.56,14.58,14.6,14.62,14.64,14.66,14.68,14.7,14.72,14.74,14.76,14.78,14.8,14.82,14.84,14.86,14.88,14.9,14.92,14.94,14.96,14.98,15,15.02,15.04,15.06,15.08,15.1,15.12,15.14,15.16,15.18,15.2,15.22,15.24,15.26,15.28,15.3,15.32,15.34,15.36,15.38,15.4,15.42,15.44,15.46,15.48,15.5,15.52,15.54,15.56,15.58,15.6,15.62,15.64,15.66,15.68,15.7,15.72,15.74,15.76,15.78,15.8,15.82,15.84,15.86,15.88,15.9,15.92,15.94,15.96,15.98,16,16.02,16.04,16.06,16.08,16.1,16.12,16.14,16.16,16.18,16.2,16.22,16.24,16.26,16.28,16.3,16.32,16.34,16.36,16.38,16.4,16.42,16.44,16.46,16.48,16.5,16.52,16.54,16.56]

pdfB = [0,0.0000086,0.000034,0.0000752,0.0001314,0.0002019,0.0002859,0.0003825,0.0004912,0.0006112,0.0007417,0.0008819,0.0010312,0.0011885,0.0013528,0.001522,0.0016919,0.0018648,0.002073,0.0022901,0.0025156,0.0027483,0.0029897,0.0032437,0.0035058,0.0037755,0.0040511,0.0043362,0.0046343,0.0049403,0.0052537,0.0055736,0.0059024,0.0062405,0.0065858,0.0069381,0.0072968,0.0076633,0.0080373,0.0084178,0.0088044,0.0091962,0.0095966,0.0100037,0.0104165,0.0108348,0.0112568,0.0116892,0.0121275,0.0125711,0.0130196,0.0134716,0.0139335,0.0144004,0.0148722,0.0153484,0.0158284,0.0163168,0.0168097,0.0173072,0.0178087,0.018314,0.0188261,0.0193423,0.0198625,0.0203864,0.020914,0.0214469,0.0219835,0.0225234,0.0230663,0.0236128,0.0241649,0.0247203,0.0252787,0.0258395,0.026404,0.0269731,0.027545,0.0281195,0.0286963,0.0292763,0.0298599,0.030446,0.0310342,0.0316243,0.0322174,0.0328133,0.033411,0.0340105,0.0346108,0.0352152,0.035822,0.0364304,0.0370401,0.0376503,0.0382645,0.0388802,0.0394971,0.040115,0.040733,0.0413551,0.0419782,0.0426021,0.0432267,0.0438513,0.0444793,0.0451079,0.0457371,0.0463664,0.0469959,0.0476281,0.0482606,0.0488934,0.049526,0.0501587,0.0507935,0.0514283,0.052063,0.0526972,0.0533315,0.0539677,0.0546036,0.0552391,0.0558737,0.0565085,0.0571444,0.0577798,0.0584143,0.0590476,0.0596814,0.0603158,0.0609492,0.0615816,0.0622124,0.0628439,0.0634753,0.0641055,0.0647345,0.0653613,0.0659893,0.0666166,0.0672425,0.0678667,0.0684884,0.0691119,0.0697339,0.0703543,0.0709729,0.0715891,0.0722061,0.0728214,0.0734348,0.0740461,0.0746552,0.0752646,0.075872,0.0764774,0.0770804,0.0776811,0.0782819,0.0788806,0.079477,0.0800707,0.0806624,0.0812534,0.081842,0.0824281,0.0830112,0.0835925,0.084173,0.0847509,0.0853262,0.0858983,0.0864687,0.0870375,0.0876036,0.0881668,0.0887265,0.0892849,0.0898413,0.0903947,0.090945,0.0914914,0.0920371,0.0925804,0.0931204,0.0936572,0.0941896,0.0947221,0.0952514,0.0957774,0.0962999,0.0968183,0.0973361,0.0978506,0.0983617,0.098869,0.0993723,0.0998748,0.1003738,0.1008691,0.1013606,0.1018482,0.1023348,0.1028177,0.1032968,0.103772,0.1042435,0.1047135,0.1051797,0.105642,0.1061,0.1065549,0.1070077,0.1074566,0.1079014,0.1083418,0.1087793,0.1092142,0.1096451,0.1100717,0.1104935,0.1109133,0.1113299,0.1117425,0.1121507,0.1125537,0.1129554,0.1133534,0.1137471,0.1141364,0.1145204,0.1149034,0.1152822,0.1156566,0.1160264,0.1163912,0.1167547,0.1171139,0.1174686,0.1178185,0.1181638,0.1185073,0.1188465,0.1191811,0.1195109,0.1198362,0.1201595,0.1204784,0.1207926,0.1211019,0.1214072,0.12171,0.1220082,0.1223018,0.1225902,0.122875,0.123157,0.1234344,0.123707,0.1239743,0.1242385,0.1244995,0.1247557,0.1250072,0.125253,0.1254965,0.1257361,0.1259709,0.1262009,0.126425,0.1266475,0.1268656,0.1270788,0.1272871,0.1274896,0.1276907,0.1278869,0.1280784,0.1282649,0.128446,0.1286252,0.1287997,0.1289693,0.1291337,0.1292932,0.1294507,0.1296032,0.1297509,0.1298933,0.1300312,0.1301667,0.1302975,0.1304232,0.1305436,0.1306601,0.1307736,0.1308825,0.1309862,0.1310845,0.1311795,0.1312713,0.1313582,0.1314402,0.1315165,0.1315902,0.13166,0.131725,0.1317848,0.1318388,0.1318912,0.1319392,0.1319822,0.1320203,0.1320525,0.1320835,0.1321095,0.1321309,0.132147,0.1321578,0.1321671,0.1321716,0.1321712,0.1321658,0.1321553,0.132143,0.1321259,0.132104,0.1320769,0.1320454,0.1320118,0.1319734,0.1319302,0.1318817,0.1318295,0.1317746,0.1317152,0.1316509,0.1315812,0.1315085,0.1314327,0.1313523,0.131267,0.1311762,0.1310833,0.1309868,0.1308857,0.1307798,0.1306679,0.1305552,0.1304383,0.1303166,0.1301901,0.130058,0.1299254,0.1297882,0.1296466,0.1295001,0.1293484,0.1291959,0.1290389,0.1288774,0.128711,0.1285402,0.128368,0.1281915,0.1280105,0.1278246,0.1276349,0.1274434,0.1272476,0.1270474,0.1268422,0.126634,0.1264235,0.1262088,0.1259897,0.1257655,0.1255393,0.1253103,0.1250772,0.1248397,0.1245969,0.1243531,0.124106,0.1238548,0.1235991,0.123338,0.1230773,0.1228125,0.1225438,0.1222707,0.1219926,0.1217146,0.1214326,0.1211467,0.1208565,0.1205619,0.120267,0.1199684,0.1196658,0.1193588,0.1190483,0.1187371,0.1184222,0.1181034,0.1177802,0.1174543,0.1171272,0.1167964,0.1164619,0.1161228,0.115782,0.1154396,0.1150936,0.1147438,0.1143894,0.1140345,0.1136773,0.1133166,0.1129522,0.112583,0.1122146,0.1118432,0.1114684,0.11109,0.1107069,0.1103251,0.10994,0.1095517,0.1091597,0.1087639,0.108369,0.1079709,0.1075696,0.1071646,0.1067566,0.1063491,0.1059387,0.105525,0.1051079,0.1046884,0.104269,0.1038466,0.1034212,0.102992,0.1025618,0.1021309,0.1016973,0.1012606,0.1008202,0.1003798,0.0999384,0.0994943,0.0990473,0.0985961,0.0981466,0.0976952,0.0972412,0.0967843,0.0963235,0.0958651,0.0954043,0.094941,0.0944749,0.0940056,0.0935384,0.0930689,0.0925971,0.0921224,0.0916454,0.0911702,0.0906928,0.0902131,0.0897307,0.0892468,0.0887642,0.0882797,0.0877929,0.0873032,0.0868133,0.0863241,0.085833,0.0853398,0.0848436,0.0843484,0.0838533,0.0833563,0.0828573,0.0823551,0.0818556,0.0813553,0.0808532,0.0803493,0.079842,0.0793387,0.0788339,0.0783275,0.0778192,0.0773085,0.0768014,0.0762929,0.075783,0.0752712,0.0747579,0.0742477,0.0737363,0.0732236,0.0727089,0.0721938,0.0716814,0.0711678,0.070653,0.0701361,0.0696202,0.0691062,0.0685913,0.0680752,0.0675568,0.0670409,0.0665262,0.0660107,0.065494,0.0649747,0.06446,0.0639454,0.06343,0.0629135,0.0623946,0.0618813,0.0613675,0.0608531,0.0603376,0.0598206,0.0593089,0.0587968,0.0582841,0.0577704,0.0572564,0.0567469,0.0562373,0.0557271,0.0552159,0.0547056,0.0541992,0.0536928,0.0531859,0.0526778,0.0521722,0.0516698,0.0511672,0.0506645,0.05016,0.0496601,0.0491624,0.0486647,0.0481668,0.0476666,0.0471737,0.0466814,0.0461893,0.045697,0.0452031,0.0447166,0.0442306,0.0437448,0.0432588,0.0427724,0.0422928,0.0418137,0.041335,0.0408559,0.0403778,0.0399058,0.0394346,0.0389637,0.0384923,0.0380235,0.03756,0.0370972,0.0366349,0.0361718,0.0357133,0.035259,0.0348055,0.0343524,0.033898,0.0334509,0.0330065,0.032563,0.0321199,0.0316751,0.0312401,0.0308063,0.0303734,0.0299409,0.0295079,0.0290842,0.0286617,0.0282403,0.0278191,0.0273988,0.0269871,0.0265767,0.0261673,0.025758,0.0253513,0.0249523,0.0245546,0.024158,0.0237609,0.0233688,0.0229831,0.0225988,0.0222155,0.0218312,0.0214547,0.021083,0.0207127,0.0203432,0.0199715,0.019612,0.0192548,0.018899,0.0185439,0.0181873,0.0178437,0.0175017,0.0171611,0.0168209,0.016481,0.016153,0.0158267,0.0155017,0.0151767,0.0148542,0.0145425,0.0142324,0.0139234,0.0136139,0.0133096,0.0130146,0.0127211,0.0124286,0.0121345,0.0118495,0.0115716,0.0112952,0.0110195,0.0107406,0.0104762,0.0102158,0.0099567,0.0096978,0.0094354,0.0091913,0.0089487,0.0087073,0.0084657,0.0082228,0.007997,0.0077726,0.007549,0.0073246,0.0071017,0.0068944,0.0066883,0.0064827,0.0062753,0.0060733,0.0058846,0.0056968,0.0055091,0.0053181,0.0051379,0.0049678,0.0047985,0.0046287,0.0044533,0.0042962,0.0041447,0.0039936,0.0038415,0.003682,0.0035479,0.0034146,0.0032814,0.0031466,0.0030077,0.0028916,0.002776,0.0026602,0.0025419,0.0024237,0.0023247,0.002226,0.0021268,0.0020242,0.0019262,0.0018433,0.0017604,0.0016766,0.0015886,0.0015101,0.0014418,0.0013733,0.0013038,0.0012286,0.0011685,0.0011131,0.0010575,0.0010006,0.0009372,0.000893,0.0008488,0.0008044,0.0007585,0.0007094,0.0006748,0.0006401,0.0006051,0.0005686,0.0005317,0.0005049,0.000478,0.0004507,0.0004218,0.0003948,0.0003743,0.0003536,0.0003326,0.0003099,0.0002907,0.0002751,0.0002594,0.0002433,0.0002253,0.0002121,0.0002003,0.0001884,0.0001761,0.000162,0.0001532,0.0001444,0.0001355,0.0001262,0.0001161,0.0001096,0.000103,0.0000963,0.0000893,0.0000822,0.0000773,0.0000725,0.0000675,0.0000622,0.0000574,0.0000538,0.0000502,0.0000466,0.0000426,0.0000393,0.0000368,0.0000342,0.0000315,0.0000285,0.0000264,0.0000246,0.0000228,0.0000209,0.0000187,0.0000174,0.0000162,0.0000149,0.0000135,0.0000121,0.0000112,0.0000103,0.0000095,0.0000085,0.0000076,0.000007,0.0000064,0.0000058,0.0000052,0.0000046,0.0000042,0.0000038,0.0000034,0.000003,0.0000027,0.0000024,0.0000022,0.0000019,0.0000016,0.0000015,0.0000013,0.0000012,0.000001,0.0000009,0.0000008,0.0000007,0.0000006,0.0000005,0.0000004,0.0000004,0.0000003,0.0000003,0.0000002,0.0000002,0.0000002,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0,0,0,0,0,0,0,0,0,0,0]

energy = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,9,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12,12.1,12.2,12.3,12.4,12.5,12.6,12.7,12.8,12.9,13,13.1,13.2,13.3,13.4,13.5,13.6,13.7,13.8,13.9,14,14.1,14.2,14.3,14.4,14.5,14.6,14.7,14.8,14.9,15,15.1,15.2,15.3,15.4,15.5]
pdf = [0.000214,0.000763,0.001513,0.002507,0.003763,0.005239,0.006914,0.008772,0.010798,0.012976,0.015292,0.017735,0.020292,0.02295,0.025699,0.028528,0.031427,0.034386,0.037395,0.040447,0.043531,0.04664,0.049767,0.052903,0.056041,0.059174,0.062296,0.065401,0.068482,0.071533,0.074549,0.077526,0.080456,0.083337,0.086164,0.088931,0.091635,0.094272,0.096839,0.099331,0.101746,0.104081,0.106332,0.108497,0.110574,0.11256,0.114452,0.11625,0.117951,0.119553,0.121056,0.122457,0.123755,0.124951,0.126042,0.127028,0.127909,0.128683,0.129351,0.129914,0.130369,0.130719,0.130963,0.131101,0.131134,0.131063,0.130888,0.130611,0.130232,0.129752,0.129174,0.128497,0.127724,0.126856,0.125895,0.124843,0.123701,0.122471,0.121156,0.119758,0.118278,0.11672,0.115086,0.113378,0.111599,0.109751,0.107838,0.105862,0.103827,0.101734,0.099587,0.09739,0.095146,0.092857,0.090528,0.088161,0.085759,0.083328,0.080869,0.078387,0.075885,0.073368,0.070837,0.068298,0.065754,0.063209,0.060667,0.058131,0.055606,0.053095,0.050602,0.048131,0.045686,0.043271,0.040889,0.038545,0.036242,0.033984,0.031774,0.029616,0.027515,0.025472,0.023493,0.021579,0.019735,0.017963,0.016266,0.014647,0.01311,0.011655,0.010286,0.009005,0.007813,0.006712,0.005703,0.004787,0.003965,0.003237,0.002602,0.002058,0.001602,0.001228,0.000929,0.000694,0.000513,0.000376,0.000273,0.000196,0.00014,0.000099,0.000069,0.000047,0.000032,0.000021,0.000014]

plt.plot(energyB, pdfB, linestyle='-', color="blue", label="J. Bahcall (1996)") 
plt.plot(energy, pdf, linestyle='-', color="orange", label="W. T. Winter (2006)") 

plt.gca().minorticks_on()
plt.gca().set_xlabel('Energy (MeV)',fontproperties=prop_font,size=32, x=1, ha='right')
plt.gca().set_ylabel(r'$\mathdefault{S}_{\nu}$',fontproperties=prop_font,size=32, y=1, ha='right')
plt.gca().get_xaxis().set_tick_params(which='both',direction='in', width=1)
plt.gca().get_yaxis().set_tick_params(which='both',direction='in', width=1)
plt.gca().yaxis.set_ticks_position('both')
plt.gca().xaxis.set_ticks_position('both')
plt.gca().get_yaxis().get_offset_text().set_fontsize(32)
plt.gca().get_yaxis().get_offset_text().set_fontproperties(prop_font)
#f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
#g = lambda rrr,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
#plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2e'))
#plt.gca().yaxis.major.formatter._useMathText = True

plt.gca().yaxis.major.formatter._useMathText = True

for label in plt.gca().get_xticklabels():
  label.set_fontproperties(prop_font)
for label in plt.gca().get_yticklabels():
  label.set_fontproperties(prop_font)

handles, labels = plt.gca().get_legend_handles_labels()
# remove the errorbars
#handles = [h[0] for h in handles]

plt.gca().legend(handles, labels,loc='upper right', fancybox=False, numpoints=1,prop=prop_font,frameon=False)

plt.show()
#exit()

# Normalize energy spectrum
sumS = 0
for i in range(len(pdf)):
  sumS = sumS + pdf[i]
for i in range(len(pdf)):
  pdf[i] = pdf[i] / sumS

sumSB = 0
for i in range(len(pdfB)):
  sumSB = sumSB + pdfB[i]
for i in range(len(pdfB)):
  pdfB[i] = pdfB[i] / sumSB

integral = 0

for i in range(len(energy)):

  ######################### Scattering Cross Sections ########################

  # Variables for the minimum and maximum Kinetic Energy of the scattered electron
  # Tmin can vary depending on the minimum limit of the energy window considered
  # Tmax is calculated using elastic scattering kinetamics
  Tmin = 0.
  Tmax = 0.

  # Global variable for the energy of the neutrino considered
  Enu = 0.

  # Cross Section arrays (for the two 7Be neutrino branches)
  sigma_nue = []
  sigma_numu = []

  # Global variable for the neutrino flavour. Options = "nue", "numu"
  type = ''

  #  for e in range(len(E_7Be)):

  Enu = energy[i]

  # Calculate Tmax
  Tmax = 0
  if Enu > 0:
    Tmax = 2. * Enu * Enu / ( me + (2.*Enu) )

    # Cross Section calculation with radiative corrections
    # First calculate cross section for the electron flavour
    type = 'nue'
    crossSect = quad(xSection, Tmin, Tmax)
    sigma_nue.append( crossSect[0] )
    #print('Nu_e cross section for E = ', str(Enu), ' MeV: ', str(crossSect[0]))

    type = 'numu'
    crossSect = quad(xSection, Tmin, Tmax)
    sigma_numu.append( crossSect[0] )
    #print('Nu_mu cross section for E = ', str(Enu), ' MeV: ', str(crossSect[0]))


    ################## Electron Neutrino Survival Probability ##################

    # Averaged Vacuum Survival Probability
    #Pee = (1 - sin2_13)**2 * (1-2*sin2_12*(1-sin2_12)) + sin2_13**2
    #print(Pee)

    # Averaged MSW Survival Probability (Adiabatic Approximation)
    # Take these values from PSelmaa
    Pee = 0.548839 -0.0230718*Enu -0.00710867*Enu*Enu + 0.0017295*Enu*Enu*Enu -0.000156206 *Enu*Enu*Enu*Enu + 6.60749e-06*Enu*Enu*Enu*Enu*Enu -1.09067e-07*Enu*Enu*Enu*Enu*Enu*Enu

    #Pee = 0.53908  -0.0221313  *Enu -0.00575418*Enu*Enu + 0.00143238*Enu*Enu*Enu -0.000129374*Enu*Enu*Enu*Enu + 5.46085e-06*Enu*Enu*Enu*Enu*Enu -8.99592e-08*Enu*Enu*Enu*Enu*Enu*Enu #Pee_+

    #Pee =  0.559208   -0.0237378*Enu -0.00872501*Enu*Enu + 0.00208222*Enu*Enu*Enu -0.00018829*Enu*Enu*Enu*Enu + 7.98985e-06*Enu*Enu*Enu*Enu*Enu -1.3228e-07*Enu*Enu*Enu*Enu*Enu*Enu  #Pee_-

    #print('Electron neutrino survival probability = ',Pee)

    integral += pdf[i] * ( Pee*sigma_nue[0] + (1-Pee)*sigma_numu[0] ) 


############################################################################
############################################################################
####
####  Calculation of the SNO+ Expected 7Be Solar Neutrino Rate
####
############################################################################
############################################################################
R = flux_7Be * n_target * integral

print(R) # events/s

R = R * 60 # events/min
print(R)

R = R * 60 # events/h
print(R)

R = R * 24 # events/d
print(R)

R = R * 365 # events/yr
print(R)

# Print the electron neutrino and muon neutrino contributions separately
#R_nue = ( flux_7Be * BR_7Be[1] * Pee[1]*sigma_nue[1] * n_target) + ( flux_7Be * BR_7Be[0] * Pee[0]*sigma_nue[0] * n_target)
#R_nue = R_nue * 60 * 60 * 24 * 365 * 7.8
#print(R_nue,' electron neutrinos /yr within the AV')

#R_numu = ( flux_7Be * BR_7Be[1] * (1-Pee[1])*sigma_numu[1] * n_target) + ( flux_7Be * BR_7Be[0] * (1-Pee[0])*sigma_numu[0] * n_target)
#R_numu = R_numu * 60 * 60 * 24 * 365 * 7.8
#print(R_numu,' muon neutrinos /yr within the AV')


#exit()


############################################################################
############################################################################
####
####  Calculation of the SNO+ Measured 8B Solar Neutrino Flux
####
############################################################################
############################################################################

frac = [0.147994 ,0.209213,0.277836]

rrr = [4.5,5,5.5]

nfit = [26.4149,30.8787,43.9561]
nfiterrp = [0.28483111,0.264076074,0.220220562] #%
nfiterrn = [0.243795301,0.229732324,0.198173972]

print('\n')
print(frac)
print('\n')
print(integral)
print('\n')
print(n_target)

print('\n Flux:')

result = []
resulterrp = []
resulterrn = []
for i in range(len(rrr)):

  result.append(nfit[i] / frac[i] / 2209.27 /60.0 /60.0 / integral / n_target)
  resulterrp.append((nfit[i]*(1.0+nfiterrp[i])) / frac[i] / 2209.27 /60.0 /60.0 / integral / n_target)
  resulterrn.append((nfit[i]*(1.0-nfiterrn[i])) / frac[i] / 2209.27 /60.0 /60.0 / integral / n_target)

print(result)

print(resulterrp)

print(resulterrn)



asymmetric_error = [resulterrn, resulterrp]
rrrr = [4.4,5,5.6]
horiz_line_data = np.array([5.66e6 for i in range(len(rrrr))])
plt.plot(rrrr, horiz_line_data, linestyle='--', color="grey", label="Observed flux, from Vitagliano et al. (2020)",linewidth=3) 

plt.errorbar(rrr, result, yerr=asymmetric_error, fmt='o', label="Fitted flux, this analysis", color="red", markersize=10)

handles, labels = plt.gca().get_legend_handles_labels()
# remove the errorbars
#handles = [h[0] for h in handles]
lgnd = plt.legend(handles,labels,loc='upper right', fancybox=False, numpoints=1, prop=prop_font,frameon=False)

plt.gca().minorticks_on()
plt.gca().set_xlabel("Fiducial volume radius (m)",fontproperties=prop_font,size=32, x=1, ha='right')
plt.gca().set_ylabel("$^{\mathdefault{8}}$B solar neutrino flux (cm$^{-\mathdefault{2}}$ s$^{-\mathdefault{1}}$)",fontproperties=prop_font,size=32, y=1, ha='right')
plt.gca().get_xaxis().set_tick_params(which='both',direction='in', width=1)
plt.gca().get_yaxis().set_tick_params(which='both',direction='in', width=1)
plt.gca().yaxis.set_ticks_position('both')
plt.gca().xaxis.set_ticks_position('both')
plt.gca().get_yaxis().get_offset_text().set_fontsize(32)
plt.gca().get_yaxis().get_offset_text().set_fontproperties(prop_font)
#f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
#g = lambda rrr,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
#plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2e'))
plt.gca().yaxis.major.formatter._useMathText = True

for label in plt.gca().get_xticklabels():
  label.set_fontproperties(prop_font)
for label in plt.gca().get_yticklabels():
  label.set_fontproperties(prop_font)

plt.xlim([4.4, 5.6])
plt.ylim([6e5, 1.8e7])

plt.show()