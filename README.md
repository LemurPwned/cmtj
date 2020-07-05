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
Please look into the [`pymtj/automation/workers.py`](pymtj/automation/workers.py) for the advanced examples of usage 
and [`test_run.py`](test_run.py) for an example on how to run test automated tasks

Please fill the issues or pull requests in case you require any additional functionality 
or encouter any problems with the software
#### Warning

Several functions allow you to set an external
Python function that updates the state of the junction
such as `.set_global_coupling_function` or `set_junction_global_external_field`  
Be very careful because setting the improper values might cause the prolonged computation times!
This library is optimized with Cython to serve the Python's `float` object which is 24 bytes. What it means for computation is that when you pass it a 32 byte object such as `numpy.float64` or even `numpy.float32`, the calculation will take significantly much more!
Take a note that for accurate results in macromagnetics you do not need such a large precision, thus you might be well off with this precision loss, expecially for the external field values.

### Coupling 

Setting the coupling value in the Junction `couplings` is defaulted to `[[2], [1]]`. Note, that the coupling is always respective to the layer ids. Layer ids start from 1 to N, where N is the last layer. Hence, `[[2], [1]]` means that the first layer is coupled with `2` layer and second layer is coupled with `1` layer.
Otherwise it can be read:
```python
[    
    [2, 3, 4] # layer 1 couples with 2, 3, 4
    [1, 3] # layer 2 couples with 1, 3,
    ... # etc. (must be provided for all layrs)
]
```

## **Calculations**

#### Conventional markings   
1. $\vec{m}_i$      the magnetisation vector of the $i$th layer. If not given, means the *current* layer 
2. $M_s$            is the magnetisation saturation 
3. $\vec{H}$        is the vector field  
#### Demagnetisation tensor 
$$\vec{H_i} = -NM_s^{j}\vec{m_i}$$
where the $M_s^{j}$ means the saturation of the *other*, *non-local* layer $j$ and $N$ is the demagnetisation tensor dependent on the elative shape of the layers $i$ and $j$.

#### Interlayer Exchange Coupling interaction (IEC)

$$\vec{H_i} = \frac{J}{\mu M_s d} (\vec{m_j} - \vec{m_i})$$
where the $J$ is the coupling constant between the layers $i$ and $j$, and $d$ is the thickness.