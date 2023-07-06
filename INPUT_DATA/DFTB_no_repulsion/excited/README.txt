#
#  Electronic Supplementary Information (ESI) parameter set for
#
#  "Data-driven approach for benchmarking DFTB-approximate excited state methods",
#   Andrés I. Bertoni and Cristián G. Sánchez*
#   *Corresponding author: csanchez@mendoza-conicet.gob.ar
#
#  This is a custom, proof-of-concept parameter set that emulates an extension of the
#  3OB parameter set [Gaus2013] to include minimal polarization on H atoms only.
#  (H now has 'p' as the maximum angular momenta.)
#
#  The repulsion splines at the end of each Slater-Koster (SK) parameter file
#  were taken directly from the original SK-files of the 3OB-3-1 parameter set [Gaus2013].
#
#  IMPORTANT: This custom parameter set should be used for excited state calculations
#  with fixed nuclei only, as we have not performed a re-parameterization of the repulsion splines.
#
#  [Gaus2013] : M. Gaus, A. Goez and M. Elstner,
#  Journal of Chemical Theory and Computation 9, 338-354 (2013).
#
#  See the rest of the ESI for the theoretical and technical details of these extensions.
#

https://pubs.rsc.org/en/content/articlelanding/2023/CP/D2CP04979A