# FFNN
A feedforward neural network (FFNN) Machine. This machine is
 constructed by providing a sequence of layers from the ``layer``
 class. Each layer implements a transformation such that the
 information is transformed sequentially as it moves from the input
 nodes through the hidden layers and to the output nodes.

## Class Constructor
Constructs a new ``FFNN`` machine:

|Argument|         Type         |            Description             |
|--------|----------------------|------------------------------------|
|hilbert |netket.hilbert.Hilbert|Hilbert space object for the system.|
|layers  |tuple                 |Tuple of layers.                    |

### Examples
A ``FFNN`` machine with 2 layers.
for a one-dimensional L=20 spin-half system:

```python
>>> from netket.layer import SumOutput
>>> from netket.layer import FullyConnected
>>> from netket.layer import Lncosh
>>> from netket.hilbert import Spin
>>> from netket.graph import Hypercube
>>> from netket.machine import FFNN
>>> g = Hypercube(length=20, n_dim=1)
>>> hi = Spin(s=0.5, total_sz=0, graph=g)
>>> layers = (FullyConnected(input_size=20,output_size=20,use_bias=True),Lncosh(input_size=20),SumOutput(input_size=20))
>>> ma = FFNN(hi, layers)
>>> print(ma.n_par)
420

```



## Class Methods 
### der_log
Member function to obtain the derivatives of log value of
machine given an input wrt the machine's parameters.

|Argument|            Type            |      Description       |
|--------|----------------------------|------------------------|
|v       |numpy.ndarray[float64[m, 1]]|Input vector to machine.|

### init_random_parameters
Member function to initialise machine parameters.

|Argument|  Type   |                               Description                                |
|--------|---------|--------------------------------------------------------------------------|
|seed    |int=1234 |The random number generator seed.                                         |
|sigma   |float=0.1|Standard deviation of normal distribution from which parameters are drawn.|

### load
Member function to load machine parameters from a json file.

|Argument|Type|             Description             |
|--------|----|-------------------------------------|
|filename|str |name of file to load parameters from.|

### log_val
Member function to obtain log value of machine given an input
vector.

|Argument|            Type            |      Description       |
|--------|----------------------------|------------------------|
|v       |numpy.ndarray[float64[m, 1]]|Input vector to machine.|

### log_val_diff
Member function to obtain difference in log value of machine
given an input and a change to the input.

|Argument|            Type            |                                 Description                                 |
|--------|----------------------------|-----------------------------------------------------------------------------|
|v       |numpy.ndarray[float64[m, 1]]|Input vector to machine.                                                     |
|tochange|List[List[int]]             |list containing the indices of the input to be changed                       |
|newconf |List[List[float]]           |list containing the new (changed) values at the indices specified in tochange|

### save
Member function to save the machine parameters.

|Argument|Type|            Description            |
|--------|----|-----------------------------------|
|filename|str |name of file to save parameters to.|

## Properties

| Property |         Type         |                                                   Description                                                    |
|----------|----------------------|------------------------------------------------------------------------------------------------------------------|
|hilbert   |netket.hilbert.Hilbert| The hilbert space object of the system.                                                                          |
|n_par     |int                   | The number of parameters in the machine.                                                                         |
|n_visible |int                   | The number of inputs into the machine aka visible units in             the case of Restricted Boltzmann Machines.|
|parameters|list                  | List containing the parameters within the layer.             Read and write                                      |
