local impl = setmetatable({}, {__index=_G})
local envmeta = {
	__index = impl,
	__newindex = function(t,k,v)
		if impl[k] then error("tried to overwrite constant") end
		print(k,v)
		impl[k] = v
	end,
}
local env = setmetatable({}, envmeta)
if setfenv then setfenv(1, env) else _ENV = env end

--[[
units
once again
using Planck units

space vs time:
speed of light: c = 1 = ~c~ m/s, for ~c~ = 2.99792458e+8
so ~c~ m = 1 s, 1/~c~ s = 1 m, m/s = 1/~c~
--]]
c = 2.99792458e+8
s_in_m = 1 / c

-- why not just change the lenght of a meter to 1/3e+8 ?
-- here is how much it would shrink by ...
--print((1/c - 1/3e+8)*c*1000)	-- difference in 1 m and 1 adj m, in terms of mm
-- 0.69180666666663 mm 

--[[
space + time vs mass:
G = 1 = ~G~ m^3 / (kg s^2), for ~G~ = 6.6740831e-11
so 1 = ~G~ m^3 / (kg s^2)
1 kg = ~G~ (m^2 / s^2) m
1 kg = ~G~/~c~^2 m
--]]
G = 6.6740831e-11
kg_in_m = G / c^2

--[[
space + time + mass vs charge:
ke = 1 = ~ke~ kg m^3 / (s^2 C^2), for ~ke~ = 8.9875517873681764e+9
so 1 = ~ke~ (~G~/~c~^2 m) m (1/~c~)^2 / C^2
C^2 = ~ke~ ~G~ / ~c~^4 m^2
1 C = sqrt(~ke~ ~G~) / ~c~^2 m
--]]
ke = 8.9875517873681764e+9
C_in_m = math.sqrt(ke * G) / c^2	-- m
N_in_m = kg_in_m / s_in_m^2	-- m^0
V_in_m = N_in_m / C_in_m	-- m^0
Ohm_in_m = kg_in_m / (s_in_m * C_in_m^2)	-- m^0

--[[
ke = 1 = 1 / (4 pi eps0) <=> eps0 = 1 / (4 pi)
c^2 = 1 / (mu0 eps0) <=> mu0 = 1/eps0 = 4 pi
--]]
mu0 = 4 * math.pi	-- m^0
eps0 = 1 / mu0	-- m^0

--[[
1 C = ~e~ e, for e = unit of charge for an electron, and ~e~ = 6.2415093414e+18
so 1 e = 1/~e~ C = sqrt(~ke~ ~G~) / (~c~^2 ~e~) m
notice that the classical electron radius is re = ke e^2 / (me c^2) = 2.817940322719e-15 m
--]]
e = 6.2415093414e+18
e_in_m = C_in_m / e

h_in_m = 1.61622938e-35
-- note the ratio of the elementary charge-converted-to-meters and Planck length 
-- is the same as the ratio of the elementary charge and Planck charge
alpha = (e_in_m/h_in_m)^2

in_in_m = .0254 

-- source: http://hyperphysics.phy-astr.gsu.edu/hbase/electric/resis.html
wire_resistivities = table{	-- at 20' Celsius, in ohm m
	aluminum = 2.65e-8,
	copper = 1.724e-8,
	iron = 9.71e-8,
	nichrome = 1e-6,
	gold = 2.24e-8,
	silver = 1.59e-8,
	platinum = 1.06e-7,
	tungsten = 5.65e-8,
}:map(function(v) return v * Ohm_in_m end)	-- ohm m => m^0
wire_diameters = table{	-- starts in inches
	electrical_range = .1019,
	household_circuit = .0808,
	switch_leads = .0640,
}:map(function(v) return v * in_in_m end)	-- in => m 

return env
