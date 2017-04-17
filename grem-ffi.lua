#!/usr/bin/env luajit
require 'ext'
require 'matrix'.__tostring = tolua	-- so when matrix_ffi defaults to matrix_lua, it will refer to tolua
local matrix = require 'matrix.ffi'
local ns = matrix{8,8,8}
print('ns',ns)
local max = matrix{1,1,1}
local min = -max
local dxs = (max - min):ediv(ns)
print('min',min)
print('max',max)
print('dxs',dxs)

local xs = matrix.lambda(ns, function(...)
	return (matrix{...} - .5):emul(dxs) + min
end)
print('xs',xs)
local eps3 = matrix{3,3,3}:lambda(function(i,j,k)
	return (i%3+1==j and j%3+1==k) and 1 or
		((k%3+1==j and j%3+1==i) and -1 or 0)
end)
-- x cross y = eps3 * y * x

