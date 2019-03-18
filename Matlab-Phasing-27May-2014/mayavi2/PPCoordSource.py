
from enthought.mayavi.sources.vtk_data_source import *
from enthought.traits.api import Instance, Array, Trait, Str, Bool, Enum, Dict
from enthought.traits.ui.api import View, Group, Item, EnumEditor, TextEditor, InstanceEditor
from enthought.mayavi.core.engine import Engine
from enthought.traits.ui.menu import Action,OKButton, CancelButton, ApplyButton
from enthought.tvtk.api import tvtk
from enthought.tvtk import array_handler
import numpy
import os
import CoordSystem3 as cs
import sys

class PPCoordSource(VTKDataSource):

    file_name=Str

    component=Enum('Amp','Phase', 'Real', 'Imag')
    space=Enum("Real", "Recip")

    transpose_input_array=Bool(True)

    coords=Instance(cs.CoordSystem)
    sg=Instance(tvtk.StructuredGrid)
    im=Instance(tvtk.ImageData)
    scalar_data=Array
    coordarray=Array
    dataarray=Array

    engine=Instance(Engine)

    set_component=Str
    scalar_name=Str

    view = View(Group(
		Item(name='transpose_input_array'),
	      	Item(name='file_name', 
			editor=TextEditor(auto_set=False,\
			enter_set=True)),
	      	Item(name='component', style='custom', editor=EnumEditor(values=component)),
#              	Item(name='scalar_name'),
#              	Item(name='vector_name'),
#             	Item(name='spacing'),
#             	Item(name='origin'),
#             	Item(name='update_image_data', show_label=False),
		label='Data'
		),

		Group(
		Item(name='coords', style='custom', show_label=False), 
		label='Coordinate System'),
#		kind='modal',
		buttons=['Ok', 'Cancel'])

    def __init__(self, **kw_args):
	super(PPCoordSource, self).__init__(**kw_args)
	fn=kw_args.pop('file_name', None)

        dir=os.path.dirname(self.file_name)
        phparams=os.path.join(dir, "phasingparams.py")

	engine=kw_args.pop('engine', None)

	self.engine.add_trait('dataarrays', Dict)
	self.engine.add_trait('coords', Dict)

	if self.engine.coords.has_key(phparams):
	   self.coords=self.engine.coords[phparams]
	else:
	   self.engine.coords[phparams]=cs.CoordSystem()
	   self.coords=self.engine.coords[phparams]
	   self.coords.exec_param_file(phparams)
	self.coords.on_trait_change(self._coords_changed, 'T')

	if fn is not None:
		self.file_name= fn
		self._open(self.file_name)

	dims=self.dataarray.shape
	print "init coords update start"
        self.coords.UpdateCoordSystem(dims)
	print "init coords update end"

	self.sg=tvtk.StructuredGrid()
	self.im=tvtk.ImageData()
	self.component="Real"
	self._component_changed('Real')		#can change this for a different default

    def _coords_changed(self, info):
	dims=self.scalar_data.shape
	self.coords.UpdateCoordSystem(dims)
	self.set_data()

    def _component_changed(self, info):
	if info=="Amp":
		self.scalar_data=numpy.abs(self.dataarray)
	if info=="Phase":
		self.scalar_data=numpy.angle(self.dataarray)
	if info=="Real":
		self.scalar_data=self.dataarray.real
	if info=="Imag":
		self.scalar_data=self.self.dataarray.imag
	self.set_data()

    def _file_name_changed(self, info):
	print self.file_name

    def _transpose_input_array_changed(self, info):
	self.set_data()

    def set_data(self):

	if not self.scalar_data.any():
		return 

	print "MAKE SGRID"
	dims=list(self.scalar_data.shape)
        self.coords.UpdateCoordSystem(dims)
	sg=self.sg
	sg.points = self.coords.coords
	sg.point_data.scalars=self.scalar_data.ravel()
	# The transpose is not needed for ScalarGrid
	if self.transpose_input_array:
	  sg.point_data.scalars=self.scalar_data.ravel()
	else:
	  sg.point_data.scalars=numpy.ravel(numpy.transpose(self.scalar_data))
	sg.point_data.scalars.name=self.component
	sg.dimensions=(dims[2], dims[1], dims[0])
	sg.extent=0, dims[2]-1, 0, dims[1]-1, 0, dims[0]-1
	sg.update_extent=0, dims[2]-1, 0, dims[1]-1, 0, dims[0]-1
	#sg.dimensions=self.scalar_data.shape
	self.data=sg
	self._update_data()
	self.update()

    
if __name__=='__main__':
	source=Sp4ArrayFileSource(file_name='/Users/rharder/PhasingProjects/pp310/22/Sequence1dist0070.sp4')
	source.configure_traits()
