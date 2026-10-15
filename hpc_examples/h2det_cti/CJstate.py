import math
import cantera as ct
from sdtoolbox.postshock import CJspeed
from sdtoolbox.postshock import PostShock_eq, PostShock_fr
from sdtoolbox.thermo import soundspeed_eq, soundspeed_fr
from optparse import OptionParser
parser = OptionParser()
OneAtm = 101300	# Pa
GasConstant = 8314
parser.add_option("--pg", type="float", dest="pg", default=101300)
parser.add_option("--te", type="float", dest="te", default=300.0)
#parser.add_option("--x", type="string", dest="comp", default='CH4:1 O2:2 Ar:5')
parser.add_option("--x", type="string", dest="comp", default='H2:2 O2:1 Ar:7')
parser.add_option("--mech", type="string", dest="mech", default='SanDiego.cti')
(init_state, args) = parser.parse_args()


# set initial quiescent state
P1 = init_state.pg;
T1 = init_state.te;

# compute initial sound speed
gas = ct.Solution(init_state.mech)
gas.TPX = T1, P1, init_state.comp
rho = gas.density
M = gas.mean_molecular_weight
cp0 = 0.001*gas.cp_mass
cv0 = 0.001*gas.cv_mass
gamma = cp0/cv0
c = math.sqrt(gamma*GasConstant*gas.T/M)
print(' ')
print('UNDISTURBED STATE')
print('Pressure = %.2f Pa' % (P1))
print('Temperature = %.2f K' % (T1))
print('Density = %.2f kg/m^3' % (rho))
print('Sound speed = %.2f m/s' % c)
print('gamma_fr = %.2f' % gamma)
print(' ')
gas()

# compute CJ Speed and CJ/VN states
[cj_speed,R2,plot_data] = CJspeed(P1, T1, init_state.comp, init_state.mech, fullOutput=True)
#[cj_speed,R2] = CJspeed(P1, T1, init_state.comp, init_state.mech, 0);   
gas_cj = PostShock_eq(cj_speed, P1, T1, init_state.comp, init_state.mech)
gas_vn = PostShock_fr(cj_speed, P1, T1, init_state.comp, init_state.mech)
Dcj = cj_speed/c
g_cj = gas_cj.cp_mass/gas_cj.cv_mass
g_vn = gas_vn.cp_mass/gas_vn.cv_mass
u_cj = cj_speed - (rho*cj_speed/gas_cj.density)
u_vn = cj_speed - (rho*cj_speed/gas_vn.density)
W_cj = gas_cj.mean_molecular_weight
W_vn = gas_vn.mean_molecular_weight
R_vn = gas_vn.cp_mass-gas_vn.cv_mass
R_cj = gas_cj.cp_mass-gas_cj.cv_mass

gas_vn.transport_model = 'Multi'
rho_ref = gas_vn.density
mu = gas_vn.viscosity
nu = gas_vn.viscosity/rho_ref





# Evaluate eq. gamma using eq. sound speed TODO give ref to Shepherd paper on eq. sound
c_cj_eq = soundspeed_eq(gas_cj)
g_cj_eq =  c_cj_eq**2*gas_cj.density/gas_cj.P




# get 65% state to find T and P
U_ps = 0.65*cj_speed	# shock speed
gas_ps = PostShock_fr(U_ps, P1, T1, init_state.comp, init_state.mech)
Ps = gas_ps.P/OneAtm
Ts = gas_ps.T
g_ps = gas_ps.cp_mass/gas_ps.cv_mass





print('CJ STATE')
print('CJ Speed = %.2f m/s ( Dcj = %.2f )' % (cj_speed,Dcj))
print('CJ Pressure = %.2f Pa' % (gas_cj.P))
print('CJ Temperature = %.2f K' % (gas_cj.T))
print('CJ Density = %.2f kg/m^3' % (gas_cj.density))
print('CJ Velocity = %.2f m/s' % u_cj)
print('CJ gamma_fr = %.2f' % g_cj)
print('CJ gamma_eq = %.2f' % g_cj_eq)
print('CJ R = %.2f' % R_cj)
print('CJ W = %.2f' % W_cj)
print('CJ Cp = %.2f' % gas_cj.cp_mass)
print(' ')
print('VN STATE')
print('VN Pressure = %.2f Pa' % (gas_vn.P))
print('VN Temperature = %.2f K' % (gas_vn.T))
print('VN Density = %.2f kg/m^3' % (gas_vn.density))
print('VN Velocity = %.2f m/s' % u_vn)
print('VN gamma_fr = %.2f' % g_vn)
print('VN R = %.2f' % R_vn)
print('VN W = %.2f' % W_vn)
print('VN Cp = %.2f' % gas_vn.cp_mass)
print('VN kin. visc. = %.3e m^2/s' % nu)
print(' ')
print('65% CJ VN STATE')
print('PS Pressure = %.2f atm' % (Ps))
print('PS Temperature = %.2f K' % (gas_ps.T))
print('PS Density = %.2f kg/m^3' % (gas_ps.density))
print('PS gamma = %.2f' % g_ps)
print(' ')


# Get heat release using VN state gamma, calc using relations from Thompson TODO clean it up
Q = c*c*(Dcj-(1/Dcj))*(Dcj-(1/Dcj))/(2*((g_vn*g_vn) - 1))
print('Heat release = %.2f J/kg' % Q)
Q = c*c*(Dcj-(1/Dcj))*(Dcj-(1/Dcj))/(2*((g_vn*g_vn) - 1))



u1 = gas_vn.int_energy_mass
uo1 = u1 - gas_vn.cv_mass*gas_vn.T
R_vn = gas_vn.cp_mass-gas_vn.cv_mass
Tvn = gas_vn.T
gas_vn.equilibrate('UV')
u2 = gas_vn.int_energy_mass
uo2 = u2 - gas_vn.cv_mass*gas_vn.T
Q = uo1 - uo2

print('CV Heat release = %.2f J/kg' % Q)
Q=Q/(R_vn*Tvn)
print('Q/RTvn (equil UV) = %.2f' % Q)




u2 = gas_cj.int_energy_mass
uo2 = u2 - gas_cj.cv_mass*gas_cj.T
Q = uo1 - uo2

print('CJ Heat release = %.2f J/kg' % Q)
Q=Q/(R_vn*Tvn)
print('Q/RTvn (CJ soln) = %.2f' % Q)






