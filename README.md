
CoMet: Combinatorial Metrics code
=================================

CoMet is an application for calculating vector similarity metrics
on large-scale parallel accelerated computing systems
to solve problems in computational genomics.
Currently the 2-way and 3-way Proportional Similarity (Czekanowski)
metrics, Custom Correlation Coefficient and DUO method are supported.
Currently the OLCF Summit system and single GPU systems are supported.
Dependencies include GCC, CUDA, MAGMA, MPI, CMake and googletest.

CoMet 1.2 New Features
----------------------

- Support for new architectures, including AMD MI100/MI250X and NVIDIA Ampere
- Support for INT4 and B1 precisions for NVIDIA Turing and Ampere architectures
- Up to 33% performance improvement for 3-way methods
- Option for lossless metrics data compression, allowing solution of much larger problems.
- Ability to set fine-grained output thresholds for different metric values.
- Option to calculate and output histograms of metric values.

Getting started
---------------

See the file Quick_Start.txt for a step-by-step guide to building and running
CoMet on the OLCF Summit system.

References
----------

W. Joubert, J. Nance, D. Weighill, D. Jacobson,
"Parallel Accelerated Vector Similarity Calculations for Genomics Applications,"
Parallel Computing, vol. 75, July 2018, pp. 130-145,
https://www.sciencedirect.com/science/article/pii/S016781911830084X,
https://arxiv.org/abs/1705.08210.

W. Joubert, J. Nance, S. Climer, D. Weighill, D. Jacobson,
"Parallel Accelerated Custom Correlation Coefficient Calculations
for Genomics Applications," Parallel Computing 84 (2019), 15-23,
https://www.sciencedirect.com/science/article/pii/S0167819118301431,
https://arxiv.org/abs/1705.08213

Wayne Joubert, Deborah Weighill, David Kainer, Sharlee Climer, Amy Justice,
Kjiersten Fagnan, Daniel Jacobson, "Attacking the Opioid Epidemic:
Determining the Epistatic and Pleiotropic Genetic Architectures
for Chronic Pain and Opioid Addiction," SC18 Gordon Bell Award paper,
https://dl.acm.org/citation.cfm?id=3291732

"GPU-enabled comparative genomics calculations on leadership-class HPC
systems," http://on-demand.gputechconf.com/gtc/2017/presentation/s7156-wayne-joubert-comparative.pdf

"CoMet: An HPC application for comparative genomics calculations,"
https://www.olcf.ornl.gov/wp-content/uploads/2017/11/2018UM-Day1-Joubert.pdf

