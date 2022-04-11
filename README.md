# 2QAN

[![PyPI version](https://badge.fury.io/py/py2QAN.svg)](https://badge.fury.io/py/py2QAN)
[![arXiv](https://img.shields.io/badge/arXiv-2108.02099-<COLOR>.svg)](https://arxiv.org/abs/2108.02099)
[![Unitary Fund](https://img.shields.io/badge/Supported%20By-UNITARY%20FUND-brightgreen.svg?style=for-the-badge)](http://unitary.fund)

2QAN is a quantum compiler for 2-local qubit Hamiltonian simulation algorithms. 2QAN uses algorithm-specific routing and scheduling techniques and can target different device topologies, different gate sets (e.g., CNOT/CX, SYC, iSWAP, sqrt_iSWAP, etc.). 

A mapping [algorithm](https://github.com/zeman412/Tabu_Search_QAP_20) based on quadratic assignment problem (QAP) is implemented in 2QAN for small benchmarks.
The QAP mapper could find good qubit initial placements for small circuits, but becomes very slow for large circuits. For fast compilation, especially for circuits with more than 40 qubits, we encourage users to use other different qubit mapping techniques such as the ones in [Qiskit](https://github.com/Qiskit/qiskit-terra/tree/main/qiskit/transpiler/passes/layout) and [tket](https://cqcl.github.io/tket/pytket/api/placement.html).

## Requirements
1. 2QAN is using Qiskit (0.36.0) circuit representation
2. Other required python packages can be found in requirement.txt

## Installation

Py2QAN can be downloaded and installed from [PyPI](https://pypi.org/project/Py2QAN/) with the command:
```
pip install py2qan
```
Note that Py2QAN requires Python 3.


## Example

Check out the examples [here](https://github.com/lllingoo/2QAN/tree/master/tests)


## Attribution

When using 2QAN for research, please cite:
```
@article{lao20212QAN,
  title={{2QAN}: A quantum compiler for 2-local qubit Hamiltonian simulation algorithms},
  author={Lao, Lingling},
  journal={arXiv preprint arXiv:2108.02099},
  year={2021}
}
```
