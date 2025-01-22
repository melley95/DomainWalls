# GRTresna

GRTresna is an open-source code for solving the constraint equations in numerical relativity.
It is developed and maintained by a collaboration of numerical relativists with a
wide range of research interests, and has already been used in many papers.

GRTresna is written entirely in C++14, using hybrid MPI/OpenMP parallelism to achieve good performance on the latest architectures.
Furthermore, it makes use of the Chombo library for adaptive mesh refinement
to allow automatic increasing of the grid resolution in regions
of arbitrary shape and topology.

Please visit www.grchombo.org for the full list of developers and their
institutions.

## Getting started
Detailed installation instructions and usage examples are available in
our [wiki](https://github.com/GRTLCollaboration/GRTresna/wiki), with the home page giving guidance on where to start.

## Contributing
We welcome feedback, bug reports, and contributions. Please consult the [wiki](https://github.com/GRTLCollaboration/GRTresna/wiki)
for our coding style and testing policy before filing a pull request.

## License
GRTresna is licensed under the BSD 3-Clause License. Please see LICENSE for details.

## Citation
If you use GRTresna as part of a paper, please cite our JOSS publication which
you can do using the following BibTeX entry:
