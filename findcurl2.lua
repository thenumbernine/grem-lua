#!/usr/bin/env luajit
require 'ext'
local matrix = require 'matrix'
local hemholtzinv = require 'matrix.hemholtzinv'
local curl = require 'matrix.curl'
matrix.__tostring = tolua
local n = matrix{8,8,8}
print('n',n)
local max = matrix{1,1,1}
local min = -max
local dx = (max - min):ediv(n)
print('min',min)
print('max',max)
print('dx',dx)
local xs = n:lambda(function(...)
	return (matrix{...} - .5):emul(dx) + min
end)

local F = n:lambda(function(...)
	local i = matrix{...}
	if i[1]==1 or i[2]==1 or i[3]==1
	or i[1]==n[1] or i[2]==n[2] or i[3]==n[3]
	then
		return matrix{0,0,0}
	end
	local xi = xs[i]
	local x,y,z = xi:unpack()
	--return matrix{-y,x,0}/xi:normSq()*.1
	return matrix{-y,x,0}
end)

local c
local F2
local function apply(F_)
	F = F_
	print('|F|',F:norm())
	c = curl(F,dx)
	print('|c|',c:norm())
	F2 = hemholtzinv{curl=c, dx=dx}
	print('|F2|',F2:norm())
	print('|F-F2|',(F-F2):norm())
end
apply(F)

local ImGuiApp = require 'imguiapp'
local View = require 'glapp.view'
local orbit = require 'glapp.orbit'
local bit = require 'bit'
local ig = require 'ffi.imgui'
local gl = require 'ffi.OpenGL'

local App = class(orbit(View.apply(ImGuiApp)))
App.viewDist = 2
function App:update()
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	App.super.update(self)

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
