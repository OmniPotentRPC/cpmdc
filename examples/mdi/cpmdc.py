"""Plugin file MDI loads for the engine name ``cpmdc``.

Point ``-plugin_path`` at this directory. A Python driver imports the
``cpmdc`` package; a C driver runs this file and calls
``MDI_Plugin_init_cpmdc``. Set ``CPMDC_PARAMS`` to the ``CPMDParams`` file.
The plugin options are ``-mdi "-role ENGINE -name CPMDC -method LINK"``.
"""

from cpmdc import MDI_Plugin_init_cpmdc

__all__ = ["MDI_Plugin_init_cpmdc"]
