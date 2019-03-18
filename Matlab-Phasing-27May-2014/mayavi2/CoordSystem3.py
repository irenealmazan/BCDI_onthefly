from enthought.traits.api import *
from enthought.traits.ui.api import *

import math as m
import numpy as nd
import os

class CoordSystem(HasTraits):

  size1 = Long(100, min=0)
  size2 = Long(100, min=0)
  size3 = Long(100, min=0)

  Wavelength = Float(-1.0)
  T=Array
  coords=Array


  Delta = Float(-1.0, auto_set=False,enter_set=True)
  Gamma = Float(-1.0, auto_set=False,enter_set=True)

  Xpixel = Float(-1, auto_set=False,enter_set=True)
  Ypixel = Float(-1, auto_set=False,enter_set=True)
  dTheta = Float(-1, min=0.0, auto_set=False,enter_set=True)

  dx = Float
  dy = Float
  dz = Float

  Arm = Float(-1, auto_set=False,enter_set=True)

  space = Enum("direct", "recip")

  view=View( 
	Item('Arm'), Item('Wavelength'),
	Item('Delta',format_str="%4.4f",springy=True),
	Item('Gamma',format_str="%4.4f"),
	Item('Xpixel',format_str="%4.4e"),Item('Ypixel',format_str="%4.4e"), 
	Item('dTheta',format_str="%4.4e",springy=True),
	Item('space'),
	resizable=True,
	buttons=['Apply','Ok', 'Cancel'] )

  def __init__(self, **traits):
	HasTraits.__init__(self, **traits)
  	oneonk=self.Wavelength/(2.0*m.pi)

#  	self.T = nd.array((1,0,0,0,1,0,0,0,1), dtype=float)
#	self.T.shape=(3,3)

  def exec_param_file(self, filename):
	if os.path.exists(filename):
		params={}
		execfile(filename, params)
		self.Xpixel = params['pixelx']
		self.Ypixel = params['pixely']
		self.Arm = params['arm']
		self.Wavelength = params['lam']
		self.Delta = params['delta']*180/m.pi
		self.Gamma = params['gam']*180/m.pi
		self.dTheta = params['dth']*180/m.pi
		self.update_transform()
	else:
		print "no params file"

  def _Delta_changed(self):
	print "delta changed"
	self.update_transform()

  def _Gamma_changed(self):
	print "gam changed"
	self.update_transform()
 
  def _space_changed(self):
	self.update_transform()

  def _Wavelength_changed(self):
	print "lam changed"
	self.update_transform()

  def _Xpixel_changed(self):
	self.update_transform() 

  def _Ypixel_changed(self):
	self.update_transform() 

  def _dTheta_changed(self):
	self.update_transform() 

  def _Arm_changed(self):
	self.update_transform() 


  def apply(self, info):
	self.update_transform()

  def update_transform(self):

	print "update trans"
	gamma=self.Gamma*m.pi/180.
	delta=self.Delta*m.pi/180.
	dth=self.dTheta*m.pi/180.

	dpx = self.Xpixel/self.Arm
	dpy = self.Ypixel/self.Arm

	dQdpx=nd.array((-m.cos(delta)*m.cos(gamma),
		0.0,
		m.sin(delta)*m.cos(gamma)))
	dQdpy=nd.array((m.sin(delta)*m.sin(gamma),
		-m.cos(gamma),
		m.cos(delta)*m.sin(gamma)))
	dQdth=nd.array((-m.cos(delta)*m.cos(gamma)+1,
		0.0,
		m.sin(delta)*m.cos(gamma)))

	Astar=(2.0*m.pi/self.Wavelength)*dpx*dQdpx
	Bstar=(2.0*m.pi/self.Wavelength)*dpy*dQdpy
	Cstar=(2.0*m.pi/self.Wavelength)*dth*dQdth

	denom = nd.dot( Astar, nd.cross(Bstar,Cstar) )
	
	A=2*m.pi*nd.cross(Bstar,Cstar)/denom
	B=2*m.pi*nd.cross(Cstar,Astar)/denom
	C=2*m.pi*nd.cross(Astar,Bstar)/denom


	if self.Delta >= 0.0 and self.Gamma >= 0.0 and self.Xpixel > 0 and \
		self.Ypixel > 0 and self.dTheta >= 0.0 and \
		self.Wavelength >= 0.0 and self.Arm >= 0.0:
	  if self.space=="direct":
		#self.T = nd.array((A,B,C)).transpose()
		self.T = nd.array((A,B,C))
	  elif self.space=="recip":
		#self.T = nd.array((Astar, Bstar, Cstar)).transpose()
		self.T = nd.array((Astar, Bstar, Cstar))
#	else:
#  	  self.T = nd.array((1,0,0,0,1,0,0,0,1), dtype=float)
#	  self.T.shape=(3,3)

	print self.T

#	print self.T[0,0], self.T[0,1], self.T[0,2]
#	print self.T[1,0], self.T[1,1], self.T[1,2]
#	print self.T[2,0], self.T[2,1], self.T[2,2]
#	print "=========================="


  def UpdateCoordSystem(self, dims):

	print "dims in getcoord", dims
	if len(dims) < 2:
		return

	if dims[0]>1:
		self.size1 = dims[0]
	if dims[1]>1:
		self.size2 = dims[1]
	if dims[2]>1:
		self.size3 = dims[2]
	else:
		self.size3 = 1
	self.dx=1./self.size1
	self.dy=1./self.size2
	self.dz=1./self.size3

	r=nd.mgrid[ (self.size1-1)*self.dx:-self.dx:-self.dx, \
		    0:self.size2*self.dy:self.dy, 0:self.size3*self.dz:self.dz]

	r.shape=3,self.size1*self.size2*self.size3
	r=r.transpose()

	self.coords = nd.dot(r, self.T)
#	print coords

	return self.coords

if __name__=="__main__":


  coord=CoordSystem()
  coord.exec_param_file("/Users/rharder/PhasingProjects/pp310/22/phasingparams.py")
  coord.GetCoordSystem( (3,3,3) )
  coord.configure_traits()

  import cstuff as c
  ccoord=c.Sp4Array()
  c.ThCoordTrans(ccoord, 3, 3, 3, .13933 ,30.5*m.pi/180, 1.0*m.pi/180, 20e-6/1.78, 40e-6/1.78, .005*m.pi/180,1.0/3,1.0/3,1.0/3, c.DIRECT)
