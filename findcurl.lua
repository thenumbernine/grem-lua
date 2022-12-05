#!/usr/bin/env luajit
require 'ext'
require 'matrix'.__tostring = tolua

local matrix = require 'matrix'
--local matrix = require 'matrix.ffi'

local n = matrix{6,6,6}
print('n',n)
local max = matrix{1,1,1}
local min = -max
local dx = (max - min):ediv(n)
print('min',min)
print('max',max)
print('dx',dx)
local xs = matrix.zeros{n[1],n[2],n[3],3}
for i in n:range() do
	xs[i] = (i - .5):emul(dx) + min
end

local F = matrix.zeros{n[1],n[2],n[3],3}
for i in n:range() do
	if i[1]==1 or i[2]==1 or i[3]==1
	or i[1]==n[1] or i[2]==n[2] or i[3]==n[3]
	then
		F[i] = matrix{0,0,0}
	end
	local xi = xs[i]
	local x,y,z = xi:unpack()
	F[i] = matrix{-y,x,0}/xi:normSq()
--[[
F = [-y,x,0]
curl F = [0, 0, 2]
--]]
	F[i] = matrix{-y,x,0}
end

local c
local F2
local function apply(F_)
	local curl = require 'matrix.curl'
	local helmholtzinv = require 'matrix.helmholtzinv'
	local startTime = os.clock()
	F = F_
	c = curl(F,dx)
	F2 = helmholtzinv{curl=c, dx=dx}
	print('took',os.clock()-startTime,'seconds')
	print('|F|',F:norm())
	print('|c|',c:norm())
	print('|F2|',F2:norm())
	print('|F-F2|',(F-F2):norm())
end
apply(F)

local bit = require 'bit'
local ig = require 'imgui'
local gl = require 'gl'

local App = require 'imguiapp.withorbit'()
App.viewDist = 2
function App:update()
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	App.super.update(self)

	gl.glColor3d(1,1,1)
	gl.glPointSize(3)
	gl.glBegin(gl.GL_POINTS)
	for i in n:range() do
		gl.glVertex3d(xs[i]:unpack())
	end
	gl.glEnd()

	gl.glBegin(gl.GL_LINES)
	for i in n:range() do
		local xi = xs[i]
		gl.glColor3d(0,1,0)
		gl.glVertex3d(xi:unpack())
		gl.glVertex3d((xi+F[i]):unpack())
		
		gl.glColor3d(0,1,1)
		gl.glVertex3d(xi:unpack())
		gl.glVertex3d((xi+F2[i]):unpack())
		
		gl.glColor3d(1,0,0)
		gl.glVertex3d(xi:unpack())
		gl.glVertex3d((xi+c[i]):unpack())
	end
	gl.glEnd()
end

function App:updateGUI()
	ig.igText('green = F')
	ig.igText('red = curl F')
	ig.igText('blue = curl^-1 curl F')

	if ig.igButton'apply' then
		apply(F2)
	end
end

App():run()
