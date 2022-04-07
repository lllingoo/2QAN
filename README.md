# 2QAN

[![Unitary Fund](https://img.shields.io/badge/Supported%20By-UNITARY%20FUND-brightgreen.svg?style=for-the-badge)](http://unitary.fund)

2QAN is a quantum compiler for 2-local qubit Hamiltonian simulation algorithms. 2QAN uses algorithm-specific routing and scheduling techniques and can target different device topologies, different gate sets. Except the QAP mapping algorithm, users can also use other different qubit mapping techniques such as the ones in [Qiskit](https://github.com/Qiskit/qiskit-terra/tree/main/qiskit/transpiler/passes/layout) and [tket](https://cqcl.github.io/tket/pytket/api/placement.html).

## Performance

## TODO
1. 2QAN is using Qiskit (0.27.0) circuit representation, it's better to use IR to avoid the pain in Qiskit updates.
2. Update tket and cirq versions, current examples with Pytket (0.11.0) and decomposition with cirq (0.11.1) 
3. Improvement in routing algorithm

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
