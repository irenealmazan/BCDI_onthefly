from enthought.mayavi.core.registry import registry
from enthought.mayavi.core.pipeline_info import PipelineInfo
from enthought.mayavi.core.metadata import ModuleMetadata

# Metadata for the new module we want to add -- notice that we use a
# factory function here for convenience, we could also use a class but
# the reasons for doing this are documented below.

#############################################################################
# READERS
#############################################################################
import spe_reader
import sp4_reader
import mat_reader

#############################################################################
#FILTERS
#############################################################################


#############################################################################
#MODULES
#############################################################################
user_outline = ModuleMetadata(
    id            = "UserOutlineModule",
    menu_name          = "&UserOutline",
    factory = 'mymod.user_outline',
    desc   = "Draw a cornered outline for given input",
    tooltip       = "Draw a cornered outline for given input",
    help       = "Draw a cornered outline for given input",
    input_info = PipelineInfo(datasets=['any'],
                              attribute_types=['any'],
                              attributes=['any'])
)

# Register the module with the mayavi registry.
registry.modules.append(user_outline)




