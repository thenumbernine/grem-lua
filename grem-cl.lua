#!/usr/bin/env luajit
require 'ext'
local ffi = require 'ffi'
local matrix = require 'matrix'
local template = require 'template'
matrix.__tostring = tolua
local ns = matrix{8,8,8}
print('ns',ns)
local max = matrix{1,1,1}
local min = -max
local dxs = (max - min):ediv(ns)
print('min',min)
print('max',max)
print('dxs',dxs)

local CLEnv = require 'cl.obj.env'
local clnumber = require 'cl.obj.number'

local function clvec(v)
	return table.map(v, function(x) 
		return tostring(clnumber(x))
	end):concat','
end

local ThisEnv = class(CLEnv)
function ThisEnv:init(...)
	ThisEnv.super.init(self, ...)
	self.code = table{self.code, template([[
#define _real3(a,b,c) (real3){.x=a, .y=b, .z=c}
constant real4 dxs = (real4)(<?=clvec(dxs)?>, 0.);
constant real4 xmin = (real4)(<?=clvec(min)?>, 0.);
constant real4 xmax = (real4)(<?=clvec(max)?>, 0.);
inline real3 real3_scale(real3 a, real s) { return _real3(a.x*s, a.y*s, a.z*s); }
inline real3 real3_cross(real3 a, real3 b) { return _real3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }

inline real eps3(int i, int j, int k) {
	if ((i+1)%3 == j && (j+1)%3 == k) return 1;
	if ((k+1)%3 == j && (j+1)%3 == i) return -1;
	return 0;
}

]], {
	clnumber = clnumber,
	clvec = clvec,
	dxs = dxs,
	min = min,
	max = max,
})}:concat'\n'
end
function ThisEnv:getTypeCode()
	return table{ThisEnv.super.getTypeCode(self), [[
typedef union {
	struct { real x, y, z; };
	struct { real s0, s1, s2; };
	real s[3];
} real3;

typedef struct {
	real s[4][4][4][4];
} Riemann_t;

typedef struct {
	real s[4][4][4];
} Conn_t;

typedef union {
	real4 m[4];
	real s[16];
} real4x4;

]]}:concat'\n'
end
local env = ThisEnv{size=ns, verbose=true}

print'xs...'
local xs = env:buffer{name='xs', type='real3'}
env:kernel{
	argsOut={xs},
	body=[[
	real4 x = ((real4)index + .5) * dxs + xmin;
	xs[index] = _real3(x.x, x.y, x.z);
]]}()

local constants = require 'constants'
local wire_radius = .5 * constants.wire_diameters.electrical_range	-- m
print('wire_radius',wire_radius)
local wire_cross_section_area = math.pi * wire_radius^2	-- m^2
print('wire_cross_section_area',wire_cross_section_area)  
local wire_length = 12 * constants.in_in_m	-- m
print('wire_length',wire_length)  
local wire_resistivity = constants.wire_resistivities.gold
print('wire_resistivity',wire_resistivity)  
local wire_resistance = wire_resistivity * wire_length / wire_cross_section_area	-- m^0
print('wire_resistance',wire_resistance)  
--local battery_voltage_in_V = 1.5
local battery_voltage_in_V = 1e+5
local battery_voltage_in_m = battery_voltage_in_V * constants.V_in_m	-- m^0
print('battery_voltage_in_m',battery_voltage_in_m)  
local wire_current = battery_voltage_in_m / wire_resistance	-- amps = C / s = m / m = m^0, likewise volts = m^0, ohms = m^0, so amps = volts / ohms = m^0
print('wire_current',wire_current)  
--local current_velocity = 1	-- doesn't matter if lambda = 0.  units of m/s = m^0
-- so ... inside the wire we know q=0 by Kirchoff's law
-- what about J = sigma E? doesn't that mean E = J rho, for rho = resistivity?
local wire_charge_density = 0	-- C / m^3 = m^-2
print('wire_charge_density',wire_charge_density)  
local wire_charge_density_per_length = wire_charge_density * wire_cross_section_area	-- m^-2 * m^2 = m^0
print('wire_charge_density_per_length',wire_charge_density_per_length)  
local wire_surface_charge_density = 0
print('wire_surface_charge_density',wire_surface_charge_density)

print'allocating cl...'
local E = env:buffer{name='E', type='real3'}
local B = env:buffer{name='B', type='real3'}
local R = env:buffer{name='R', type='Riemann_t'}
local dt_Conn = env:buffer{name='dt_Conn', type='Conn_t'}
local div_Conni = env:buffer{name='div_Conni', type='real4x4'}
local Conn = env:buffer{name='Conn', type='Conn_t'}

-- Conn^a_bt,i: [x,y,z,i] for fixed a,b
-- the calculations stored in this can be merged into the lapConnUa_bt kernel if we need the memory
local ConnUa_bt_i = env:buffer{name='ConnUa_bt_i', type='real3'}

-- Conn^a_bt,ij g^ij: [x,y,z] for fixed a,b
-- using g^ij = delta^ij for now
local lapConnUa_bt = env:buffer{name='lapConnUa_bt', type='real'}

print'building kernels...'
local init_E_B = env:kernel{
	argsOut={E,B},
	argsIn={xs},
	body=template([[
	real3 xi = xs[index];
	real x = xi.x;
	real y = xi.y;
	real z = xi.z;	// m
	real r2 = sqrt(x*x+y*y);	// m
	
	// hyperphysics:
	//local lambda = wire_charge_density_per_length / current_velocity	-- m^0 / m^0 = m^0
	//local Er = lambda / (2 * math.pi) / (constants.eps0 * r2)	-- m^0 / m = m^-1
	// https://physics.stackexchange.com/questions/291779/electric-field-outside-wire-with-stationary-current?rq=1 
	// J = sigma E
	real Ez_int = <?=wire_current * wire_resistivity?>;	// m^0 * m^0 = m^0
	real Er_int = 0.;
	real Er = <?=wire_surface_charge_density * wire_radius?> / (<?=constants.eps0?> * r2);	// m^0 * m / m = m^0
	real Ez = <?=wire_current * wire_resistivity?>;	// m^0
	E[index] = real3_scale(_real3(x/r2*Er,y/r2*Er,Ez), <?=math.sqrt(constants.eps0)?>);
	// times sqrt(eps0) for convenience of representing T_ab

	// http://www.ifi.unicamp.br/~assis/Found-Phys-V29-p729-753(1999).pdf
	real Bt_int = <?=constants.mu0 * wire_current?> * r2 / <?=2 * math.pi * wire_radius^2?>; // m^0 m / (m^2) = m^-1
	real Bt = <?=constants.mu0 * wire_current?> / (<?=2 * math.pi?> * r2);	// m^-1
	//print('Bt',Bt)
	B[index] = real3_scale(_real3(-y/r2 * Bt, x/r2 * Bt, 0), <?=1/math.sqrt(constants.mu0)?>);
	// divide sqrt(mu0) for convenience ... same deal
]], {
	wire_current = wire_current,
	wire_resistivity = wire_resistivity,
	wire_surface_charge_density = wire_surface_charge_density,
	wire_radius = wire_radius,
	constants = constants,
})}
init_E_B:compile() 

local init_R = env:kernel{
	argsOut={R},
	argsIn={E,B},
	body=[[
	real3 Ei = E[index];
	real3 Bi = B[index];
	real3 Si = real3_cross(Ei,Bi);
	global Riemann_t* Ri = R+index;
	for (int a = 0; a < 4; ++a) {
		for (int b = 0; b < 4; ++b) {
			for (int c = 0; c < 4; ++c) {
				for (int d = 0; d < 4; ++d) {
					real s = 1;
					real tmp;
					if (a>b) {
						tmp = a; a = b; b = tmp; s = -s;
					}
					if (c>d) {
						tmp = c; c = d; d = tmp; s = -s;
					}
					if (a==0 && b>0 && c==0 && d>0) {
						// R^t_itj = -E_i E_j - B_i B_j
						Ri->s[a][b][c][d] = -s * (Ei.s[b-1] * Ei.s[d-1] + Bi.s[b-1] * Bi.s[d-1]);
					} else if (a==0 && b>0 && c>0 && d>0) {
						// R^t_ijk = gamma_ij S_k - gamma_ik S_j
						Ri->s[a][b][c][d] = s * (
							(b==c ? Si.s[d-1] : 0) 
							- (b==d ? Si.s[c-1] : 0)
						);
					} else if (a>0 && b>0 && c>0 && d>0) {
						// R^i_jkl = eps^i_jm (E^m E_n + B^m B_n) eps^n_jk
						real sum = 0;
						for (int m = 0; m < 3; ++m) {
							for (int n = 0; n < 3; ++n) {
								sum += eps3(a-1,b-1,m) * (Ei.s[m] * Ei.s[n] + Bi.s[m] * Bi.s[n]) * eps3(n,c-1,d-1);
							}
						}
						Ri->s[a][b][c][d] = sum;
					} else {
						Ri->s[a][b][c][d] = 0;
					}
				}
			}
		}
	}
]]}
init_R:compile() 

local calc_ConnUa_bt_i = range(4):map(function(a)
	return range(4):map(function(b)
		local k = env:kernel{
			argsOut = {ConnUa_bt_i},
			argsIn = {Conn, dt_Conn},
			body = template([[
	const int a = <?=a-1?>;
	const int b = <?=b-1?>;
	real ConnSq = 0;
	global const Conn_t* Conni = Conn + index;
	for (int i = 0; i < 3; ++i) {
		for (int c = 0; c < 4; ++c) {
			ConnSq += Conni->s[a][c][b] * Conni->s[c][b][i+1]
					- Conni->s[a][c][i+1] * Conni->s[c][b][0];
		}
	}
]], 	{
			a = a,
			b = b,
		})}
		k:compile()
		return k
	end)
end)

local calc_lapConnUa_bt = range(4):map(function(a)
	return range(4):map(function(b)
		local k = env:kernel{
			argsOut = {lapConnUa_bt},
			argsIn = {ConnUa_bt_i},
			body = template[[
	if (i.x == 0 || i.x == size.x-1 || 
		i.y == 0 || i.y == size.y-1 ||
		i.z == 0 || i.z == size.z-1) 
	{
		lapConnUa_bt[index] = 0;
	} else {
		lapConnUa_bt[index] = 0.
<? for j=0,2 do 
?>			+ .5 * dxs.s<?=j?> * (
				ConnUa_bt_i[index+stepsize.s<?=j?> ].s<?=j?> - 
				ConnUa_bt_i[index-stepsize.s<?=j?> ].s<?=j?>)
<? end 
?>		;
	} 
]]}
		k:compile()
		return k
	end)
end)

print'executing...'
print'E and B...'
init_E_B()

print'R...'
init_R()

--dt_Conn:fill()
--div_Conni:fill()
--Conn:fill()

for iter=1,20 do
	print('iter',iter)
	for a=1,4 do
		for b=1,4 do
			-- store Conn^a_bt,ij g^ij
			print('calc_ConnUa_bt_i',a,b)
			calc_ConnUa_bt_i[a][b]()

			-- calculate Conn^a_bt
			print('calc_lapConnUa_bt',a,b)
			calc_lapConnUa_bt[a][b]()
		end
	end
end

print'done'
