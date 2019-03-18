"""A simple source that allows one to view a suitably shaped numpy
array as ImageData.  This supports both scalar and vector data.
"""
# Author: Prabhu Ramachandran <prabhu_r@users.sf.net>
# Copyright (c) 2005, Enthought, Inc.
# License: BSD Style.

# Standard library imports.
import numpy

# Enthought library imports
from enthought.traits.api import Instance, Array, Trait, Str, Bool, Button
from enthought.traits.ui.api import View, Group, Item
from enthought.tvtk.api import tvtk
from enthought.tvtk import array_handler

# Local imports
from enthought.mayavi.core.source import Source
from enthought.mayavi.core.pipeline_info import PipelineInfo

def _check_scalar_array(obj, name, value):
    """Validates a scalar array passed to the object."""
    if value is None:
        return None
    arr = numpy.asarray(value)
    assert len(arr.shape) in [2,3], "Scalar array must be 2 or 3 dimensional"
    return arr

_check_scalar_array.info = 'a 2D or 3D numpy array'



######################################################################
# 'ArraySource' class.
######################################################################
class FitsSource(Source):

    """A simple source that allows one to view a suitably shaped numpy
    array as ImageData.  This supports both scalar and vector data.
    """

    # The scalar array data we manage.
    scalar_data = Trait(None, _check_scalar_array, rich_compare=False)

    # The name of our scalar array.
    scalar_name = Str('scalar')
    
    # The spacing of the points in the array.
    spacing = Array(dtype=float, shape=(3,), value=(1.0, 1.0, 1.0),
                    desc='the spacing between points in array')

    # The origin of the points in the array.
    origin = Array(dtype=float, shape=(3,), value=(0.0, 0.0, 0.0),
                    desc='the origin of the points in array')

    # Fire an event to update the spacing and origin - this reflushes
    # the pipeline.
    update_image_data = Button('Update spacing and origin')

    # The image data stored by this instance.
    image_data = Instance(tvtk.ImageData, allow_none=False)

    # Should we transpose the input data or not.  Transposing is
    # necessary to make the numpy array compatible with the way VTK
    # needs it.  However, transposing numpy arrays makes them
    # non-contiguous where the data is copied by VTK.  Thus, when the
    # user explicitly requests that transpose_input_array is false
    # then we assume that the array has already been suitably
    # formatted by the user.
    transpose_input_array = Bool(True, desc='if input array should be transposed (if on VTK will copy the input data)')

    # Information about what this object can produce.
    output_info = PipelineInfo(datasets=['image_data'])

    # Our view.
    view = View(Group(Item(name='transpose_input_array'),
                      Item(name='scalar_name'),
                      Item(name='spacing'),
                      Item(name='origin'),
                      Item(name='update_image_data', show_label=False),
                      show_labels=True)
                )
    
    ######################################################################
    # `object` interface.
    ######################################################################
    def __init__(self, **traits):
        # Set the scalar and vector data at the end so we pop it here.
        sd = traits.pop('scalar_data', None)
        # Now set the other traits.
        super(FitsSource, self).__init__(**traits)
        # And finally set the scalar and vector data.
        if sd is not None:
            self.scalar_data = sd

        # Setup the mayavi pipeline by sticking the image data into
        # our outputs.
        self.outputs = [self.image_data]

    def __get_pure_state__(self):
        d = super(FitsSource, self).__get_pure_state__()
        d.pop('image_data', None)
        return d
    
    ######################################################################
    # ArraySource interface.
    ######################################################################
    def update(self):
        """Call this function when you change the array data
        in-place."""
        d = self.image_data
        d.modified()
        pd = d.point_data
        if self.scalar_data is not None:
            pd.scalars.modified()
        self.data_changed = True

    ######################################################################
    # Non-public interface.
    ######################################################################
    def _image_data_default(self):
        s = tuple(self.spacing)
        o = tuple(self.origin)
        return tvtk.ImageData(spacing=s, origin=o)

    def _image_data_changed(self, value):
        self.outputs = [value]

    def _update_image_data_fired(self):
        sp = tuple(self.spacing)
        o = tuple(self.origin)
        self.image_data = tvtk.ImageData(spacing=sp, origin=o)
        sd = self.scalar_data
        if sd is not None:
            self._scalar_data_changed(sd)
    
    def _scalar_data_changed(self, data):
        img_data = self.image_data
        if data is None:
            img_data.point_data.scalars = None
            self.data_changed = True
            return
        dims = list(data.shape)
        if len(dims) == 2:
            dims.append(1)
      
        img_data.origin = tuple(self.origin)
        img_data.dimensions = tuple(dims)
        img_data.extent = 0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1
        img_data.update_extent = 0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1
        if self.transpose_input_array:
            img_data.point_data.scalars = numpy.ravel(numpy.transpose(data))
        else:
            img_data.point_data.scalars = numpy.ravel(data)            
        img_data.point_data.scalars.name = self.scalar_name
        # This is very important and if not done can lead to a segfault!
        typecode = data.dtype
        img_data.scalar_type = array_handler.get_vtk_array_type(typecode)
        img_data.update() # This sets up the extents correctly.
        img_data.update_traits()

        # Now flush the mayavi pipeline.
        self.data_changed = True

    def _scalar_name_changed(self, value):
        if self.scalar_data is not None:
            self.image_data.point_data.scalars.name = value
            self.data_changed = True
            
    def _transpose_input_array_changed(self, value):
        if self.scalar_data is not None:
            self._scalar_data_changed(self.scalar_data)
