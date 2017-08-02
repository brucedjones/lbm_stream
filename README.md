# LBM_Stream

This a sample implementation of Lagrangian/Origin Streaming for the Lattice Boltzmann Method.

Using preprocessor directives, the code may be compiled using either origin streaming or the traditional two-lattice approach to streaming.

## Build instructions
Prerequisites: Cmake

```bash
git clone <github_url>
mkdir build && cd build
cmake <options> ../lbm_stream
make
./lbm_stream
```

Options:
  * -DORIGIN_STREAMING - Origin streaming implementation
  * -DTWO_LATTICE - Traditional two lattice approach
  * -DOUTPUT - Executeable will output a binary file of LX x LY floats representing x-velocity for each lattice site

Note: Do not use ORIGIN_STREAMING and TWO_LATTICE options at the same time


