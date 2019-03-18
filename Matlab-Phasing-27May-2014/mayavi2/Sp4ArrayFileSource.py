
from enthought.mayavi.sources.api import ArraySource
from enthought.traits.api import Instance, Array, Trait, Str, Bool, Enum
from enthought.traits.ui.api import View, Group, Item, EnumEditor, TextEditor
import numpy
import Sp4File as sp4

class Sp4ArrayFileSource(ArraySource):

    file_name=Str

    component=Enum('Amp','Phase', 'Real', 'Imag')
    set_component=Str

    f=Instance(sp4.Sp4File)

    view = View(Group(Item(name='transpose_input_array'),
		      Item(name='file_name', editor=TextEditor(auto_set=False,\
				  			       enter_set=True)),
		      Item(name='component', editor=EnumEditor(values=component)),
                      Item(name='scalar_name'),
                      Item(name='vector_name'),
                      Item(name='spacing'),
                      Item(name='origin'),
                      Item(name='update_image_data', show_label=False),
                      show_labels=True)
                )


    def __init__(self, **kw_args):
	super(Sp4ArrayFileSource, self).__init__(**kw_args)
	fn=kw_args.pop('file_name', None)
	if fn is not None:
		self.file_name= fn
		self._open(self.file_name)
	self.component="Amp"
	self._component_changed('Amp')

    def _open(self, fn):
	self.f=sp4.Sp4File(self.file_name)

    def _component_changed(self, info):
	if info=="Amp":
		self.scalar_data=numpy.abs(self.f.GetArray())
	if info=="Phase":
		self.scalar_data=numpy.angle(self.f.GetArray())
	self.update()

    def _file_name_changed(self, info):
	print self.file_name

    
if __name__=='__main__':
	source=Sp4ArrayFileSource(file_name='/Volumes/Disk1/PhasingProjects/MoyuGoldDifferenceMap/BD-data81.sp4')
	source.configure_traits()
