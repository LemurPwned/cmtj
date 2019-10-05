### SPINPY

A Python library for spin magnetics

### Warning

Several functions allow you to set an external
Python function that updates the state of the junction
such as `.set_global_coupling_function` or `set_junction_global_external_field`  
Be very careful because setting the improper values might cause the prolonged computation times!
This library is optimized with Cython to serve the Python's `float` object which is 24 bytes. What it means for computation is that when you pass it a 32 byte object such as `numpy.float64` or even `numpy.float32`, the calculation will take significantly much more!
Take a note that for accurate results in macromagnetics you do not need such a large precision, thus you might be well off with this precision loss, expecially for the external field values.
