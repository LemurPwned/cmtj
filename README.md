## PYMTJ

A Python library for spin magnetics in Magnetic Tunnel Junctions.

### Installation
Installation is simple and requires Python3 (preferably 3.7)
Clone this repository
```bash 
git clone https://github.com/LemurPwned/spinpy.git
```
Enter the repository and: 
```bash 
pip3 install -e .
```
This will install and the dependencies automatically and all the dependencies as well.
For the precise requirements that will be installed please see `setup.py` files and the 
requirements section.


### Using the software
Please look into the `pymtj\automation\workers.py` for the advanced examples of usage 
and `test_run.py` for an example on how to run test automated tasks

Please fill the issues or pull requests in case you require any additional functionality 
or encouter any problems with the software
#### Warning

Several functions allow you to set an external
Python function that updates the state of the junction
such as `.set_global_coupling_function` or `set_junction_global_external_field`  
Be very careful because setting the improper values might cause the prolonged computation times!
This library is optimized with Cython to serve the Python's `float` object which is 24 bytes. What it means for computation is that when you pass it a 32 byte object such as `numpy.float64` or even `numpy.float32`, the calculation will take significantly much more!
Take a note that for accurate results in macromagnetics you do not need such a large precision, thus you might be well off with this precision loss, expecially for the external field values.
