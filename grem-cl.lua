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
]]}:concat'\n'
end
local env = ThisEnv{size=ns}

local xs = env:buffer{name='xs', type='real3'}
env:kernel{
	argsOut={xs},
	body=[[
	real4 x = ((real4)index + .5) * dxs + xmin;
	xs[index] = _real3(x.x, x.y, x.z);
]]}()

local constants = require 'constants'
