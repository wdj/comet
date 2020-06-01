
CoMet: Combinatorial Metrics code
=================================

CoMet is an application for calculating vector similarity metrics
on large-scale parallel accelerated computing systems
to solve problems in computational genomics.
Currently the 2-way and 3-way Proportional Similarity (Czekanowski)
metrics, Custom Correlation Coefficient and DUO method are supported.
Currently the OLCF Summit system and single GPU systems are supported.
Dependencies include GCC, CUDA, MPI, CMake and googletest.

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
for Genomics Applications," arxiv 1705.08213 [cs], Parallel Computing, accepted,
https://arxiv.org/abs/1705.08213

Wayne Joubert, Deborah Weighill, David Kainer, Sharlee Climer, Amy Justice,
Kjiersten Fagnan, Daniel Jacobson, "Attacking the Opioid Epidemic:
Determining the Epistatic and Pleiotropic Genetic Architectures
for Chronic Pain and Opioid Addiction," SC18 Gordon Bell paper,
https://dl.acm.org/citation.cfm?id=3291732

"GPU-enabled comparative genomics calculations on leadership-class HPC systems,"
http://on-demand.gputechconf.com/gtc/2017/presentation/s7156-wayne-joubert-comparative.pdf

"CoMet: An HPC application for comparative genomics calculations,"
https://www.olcf.ornl.gov/wp-content/uploads/2017/11/2018UM-Day1-Joubert.pdf

Wayne Joubert, Deborah Weighill, David Kainer, Sharlee Climer, Amy Justice,
Kjiersten Fagnan, Daniel Jacobson, "Attacking the Opioid Epidemic:
Determining the Epistatic and Pleiotropic Genetic Architectures
for Chronic Pain and Opioid Addiction."  Supercomputing SC18, Gordon Bell Award.

