Documentation for the program SCATCI_INTEGRALS
==============================================
* 13/09/2016 ZM
* 09/11/2018 ZM

SCATCI_INTEGRALS computes the molecular integrals needed in the rest of the UKRmol+ calculation.
The target GTO basis and the molecular orbitals are loaded from a formatted file in the Molden format (this can be generated for example using Molpro or Gaussian).
The basis for the continuum can comprise GTOs, BTOs or both and is specified on input to SCATCI_INTEGRALS.

The program performs orthogonalization of the continuum basis against the target orbitals, deletes the continuum orbitals using the specified threshold and then proceeds to evaluate all 1-electron and 2-electron integrals in the basis of the target and the continuum orbitals.
The molecular integrals evaluated can be restriced to 1 particle in the continuum.

This program can be compiled in double or in quad precision. The quad precision version typically allows to retain the full continuum basis without problems with linear dependencies.

Usage notes
------------

SCATCI_INTEGRALS requires namelists

* `&target_data ... /`

* `&continuum_data ... /`

* `&process_control ... /`

and takes the entirely optional namelist

* `&pco_data ... /`


The format of the namelist variables described below is the following:

NAME TYPE DIMENSION (default value)

Namelist `target_data`
--------------------

* **molden_file C 1 ('')**

   Path to the Molden file containing the molecular geometry, GTO basis and the molecular orbitals.

* **nob I 8 (0)**

   Number of molecular orbitals per symmetry to load from the Molden file.

* **no_sym_op I 1 (0)**

   Number of symmetry operations to use.

* **sym_op C 3 (' ')**

   Symmetry operation flags. Exactly no_sym_op flags must be given. These can be one of: X, Y, Z, XY, YZ, XZ, XYZ and determine which axes change sign under the symmetry operation.
   The symmetry operations applied must be consistent with the symmetry of the molecular orbitals saved on the Molden file.

* **a R 1 (-1.0)**

   R-matrix radius. A value .le. 0 is equivalent to infinity.

* **select_orbitals_by I 1 (1)**

   If set to 1 then the array nob selects orbitals according to their index from the Molden file (default). If set to 2 then nob(i)=x selects the x lowest-lying orbitals from symmetry i.
   The orbital energy is taken from the Ene flag of the orbitals on the Molden file.

* **alpha_or_beta I 1 (1)**

   Select only orbitals of the chosen spin from the Molden file. The options are: 0 = allow both spins, 1 = use alpha orbitals only,
   2 = use beta orbitals only.

Namelist `continuum_data`
-----------------------

* **del_thrs R 8 (-1.0)**

   Threshold for deletion of the continuum orbitals from each symmetry after symmetric orthogonalization.

* **min_l I 1 (-1)**

   The lowest GTO continuum partial wave to include in the calculation.

* **max_l I 1 (-1)**

   The highest GTO continuum partial wave to include in the calculation.

* **exponents R (1:20,0:15) (0.0)**

   Exponents of the continuum GTOs for each partial waves. Exponents for partial wave angular momenta L = [min_l,max_l] will be included in the calculation. Exponents for partial waves outside of this range are ignored if specified. E.g. GTO continuum basis for partial waves L = [0,6] is specified as:

   ```fortran
   min_l = 0, max_l = 6,
   exponents(:,0) = 33.559323, 2.576968, 0.276187, 0.184632, 0.126330, 0.086941, 0.059603, 0.040429, 0.026977, 0.017594, 0.011082,
   exponents(:,1) =  1.546842, 0.367220, 0.254952, 0.182860, 0.133060, 0.097244, 0.070868, 0.051199, 0.036466, 0.025435, 0.017157,
   exponents(:,2) =  0.437039, 0.303223, 0.217211, 0.158560, 0.116878, 0.086332, 0.063474, 0.046155, 0.032943, 0.022774,
   exponents(:,3) =  0.409883, 0.291455, 0.213760, 0.159587, 0.120199, 0.090678, 0.068084, 0.050566, 0.036876, 0.026055,
   exponents(:,4) =  0.340244, 0.243185, 0.179271, 0.134479, 0.101646, 0.076790, 0.057555, 0.042443, 0.030365,
   exponents(:,5) =  0.225709, 0.167874, 0.128505, 0.099316, 0.076719, 0.058780, 0.044303, 0.032392,
   exponents(:,6) =  0.118188, 0.085985, 0.062824, 0.045463, 0.032349, 0.022480, 0.015073,
   ```

   If continuum basis is to be included then both min_l and max_l must be specified. If neither of min_l,max_l are set to >= 0 then the GTO continuum basis is not included.

* **min_bspline_l I 1 (-1)**

   The lowest BTO continuum partial wave to include in the calculation.

* **max_bspline_l I 1 (-1)**

   The highest BTO continuum partial wave to include in the calculation.

* **bspline_grid_start R 1 (-1.0)**

   Radial distance from the center of coordinates where the B-spline basis starts (it always ends on r=a, i.e. on the R-matrix radius).

* **bspline_order I 1 (-1)**

   Order of the radial B-splines to include in the calculation (typically in the range 8 to 11).

* **no_bsplines I 1 (-1)**

   The size of the B-spline basis, i.e. number of radial B-splines to generate.

* **bspline_indices I (2,0:15) (-1)**

   For each BTO angular momentum specify the first and the last B-spline to include. E.g. to specify that for L=2 BTOs the radial B-splines with indices 2 to 20 should be included we set bspline_indices as:

   ```fortran
   bspline_indices(1,2) = 3,
   bspline_indices(2,2) = 20,
   ```

   If the B-spline basis starts in the center of coordinates then for partial wave with angular momentum L B-splines with indices < (L+1) should not be included in the basis due to boundary condition 
   requirements (r -> 0: psi(r) ~ r**(L+1)). If the BTO basis starts at r > 0.0 then the first two radial B-splines must not be included (these have a non-zero first derivative at the starting point).
   
* **run_free_scattering L 1 (.false.)**

   If a continuum basis is included in the calculation this flag allows to run the free-potential scattering calculation (useful for testing quality of the continuum basis).

* **min_energy R 1 (0.0)**

   Starting energy of the continuum electron for the free-potential scattering calculation.

* **max_energy R 1 (1.0)**

   Final energy of the continuum electron for the free-potential scattering calculation.

* **nE I 1 (-1)**

   Number of energy points in the range [min_energy,max_energy] for the free-potential scattering calculation.

Namelist `process_control`
------------------------

* **basis_input C 1 ('')**

   Path to a GBTOlib datafile (e.g. earlier generated 'moints' file) with a basis set to use. The basis setup in other namelists needs to be identical, but the basis will be read from the file
   rather than freshly generated. This is useful when sign-compatible basis sets are needed for calculation of another orbital properties (e.g. damped dipoles).

* **max_ijrs_size R 1 (-1.0)**

   The maximum allowed size (in Mib) of an intermediate array used in the 2-electron integral transformation. The larger the better: influences the speed of the transformation.

* **save_ao_integrals_to_disk L 1 (.false.)**

   Flag to force saving of the atomic orbital integrals to disk. At the moment this functionality can be used but is not useful since loading of the atomic integrals and subsequent transformation has not been 
   implemented in scatci_integrals. However, this functionality exists in the integral library so the integral transformation starting from a pregenerated file with atomic integrals can be implemented as a separate program.

* **ao_integrals_file_name C 1 ('aoints')**

   Name of the file containing the atomic orbital integrals.

* **mo_integrals_file_name C 1 ('moints')**

   Name of the file containing the molecular orbitals, atomic basis and the transformed molecular orbital integrals.

* **redirect_master L 1 (.true.)**

   Redirect the text output of the master process to a text file. Text outputs of other ranks are always redirected, regardless of this flag.

* **do_two_particle_integrals C 1 (.true.)**

   Flag to select calculation of 2-electron integrals or not.

* **use_spherical_cgto_alg C 1 (.true.)**

   Flag to select spherical-coordinate-based algorithms for generation of the GTO-only integrals. The opposite setting is strongly discouraged since it selects cartesian-based algorithms
   which in the present implementation can suffer from numerical instabilities.

* **check_target_target_orbital_overlaps L 1 (.true.)**

   Flag to enable checking of the overlaps between the target orbitals against a hard-coded threshold as defined in the module const.f90.

* **check_target_continuum_orbital_overlaps L 1 (.true.)**

   Flag to enable checking of the overlaps between the target and the continuum orbitals against a hard-coded threshold as defined in the module const.f90.

* **check_continuum_continuum_orbital_overlaps L 1 (.true.)**

   Flag to enable checking of the overlaps between the continuum orbitals against a hard-coded threshold as defined in the module const.f90.

* **two_p_continuum L 1 (.false.)**

   Flag to enable calculation of integrals for 2 electrons in the continuum. For UKRmol+ we should always use .false. since we don't need integrals for 2p in the continuum: we gain a significant speed up the integral evaluation. 

* **mixed_ints_method I 1 (-1)**

   Method for calculation of the mixed BTO/GTO integrals:  
   = 1: Legendre expansion  
   = 2: Lebedev quadrature  
   = 3: Semi-analytic for 1-electron integrals and Legendre expansion for 2-electron integrals (recommended)

* **molecular_2el_algorithm I 1 (0)**

   Method for transformation of atomic two-electron integrals to molecular two-electron integrals.  
   = 0: Sparse transformation if continuum basis consists only of B-splines, dense transformation otherwise.  
   = 1: Sparse transformation.  
   = any other: Dense transformation.

* **max_l_legendre_1el I 1 (-1)**

   Maximum L in the Legendre expansion of the mixed nuclear attraction integrals (typically > 50). Applicable if mixed_ints_method = 1.

* **max_l_legendre_2el I 1 (-1)**

   Maximum L in the Legendre expansion of the mixed 2-electron integrals (typically > 50).

* **scratch_directory C 1 ('')**

   Directory where temporary files can be stored. If this path is specified then the auxiliary Y_lm functions needed for the mixed integral evaluation will be saved to this directory instead of keeping them in the memory. This makes sense to do only if the disk storage is fast (e.g. SSD). At the moment the scratch_directory is not used in any other part of the code.

* **delta_r1 R 1 (0.25)**

   The length (in Bohr) of the primitive quadrature interval for evaluation of the mixed BTO/GTO integrals. The smaller the interval the larger the compute time.

* **calc_radial_densities L 1 (.false.)**

   If set to .true. then the radial charge densities for all orbitals in the basis will be calculated.

* **print_1el_ints L 1 (.false.)**

   If set to .true. then the 1-electron atomic integrals will be printed.

* **print_2el_ints L 1 (.false.)**

   If set to .true. then the 2-electron atomic and molecular integrals will be printed.

* **preorthogonalize_continuum L 1 (.false.)**

   Orthogonalize the CGTO and BTO continuum separately first.

* **ortho_continuum_against_all_tgt_orbs L 1 (.false.)**

   If set, the continuum will be orthogonalized against the full set of target orbitals. Otherwise the continuum will be
   orthogonalized as usual against the number of `nob(:)` orbitals from each symmetry.

* **qmoln L 1 (.false.)**

   Internal cooperation flag used by the Quantemol EC software.

* **verbosity I 1 (1)**

   Text output verbosity level. Meaningful values are 0 for no text output at all, 1 (default) for output of the main information, 2 for extended information and 3 for debugging.

* **dipole_damp_factor R 1 (0.0)**

   When non-zero the dipole properties are calculated with the exponentially damped radial part of the dipole operator: `r * exp(-dipole_damp_factor * r)`.
   The value of the damping factor is saved on the moints file in the structure integral_options.


Sample input file
-----------------

Input file used for a large photoionization calculation for the CO2 molecule in GTO basis.
This calculation was performed in quad precision (note the low value of `del_thrs`).

```fortran
&target_data
 a = 10.0,
 no_sym_op = 3,
 sym_op = 'X', 'Y', 'Z',
 molden_file =  './co2.molden'
 nob = 11, 7, 6, 2, 9, 5, 4, 2,
/
&continuum_data
  run_free_scattering = .true.,  min_energy = 0.0, max_energy = 3.0, nE = 300,
  del_thrs = 1.0D-14,1.0D-14,1.0D-14,1.0D-14,1.0D-14,1.0D-14,1.0D-14,1.0D-14,
  min_l = 0, max_l = 7,
  exponents(:,0) = 33.559323, 2.576968, 0.276187, 0.184632, 0.126330, 0.086941, 0.059603, 0.040429, 0.026977, 0.017594, 0.011082,
  exponents(:,1) =  1.546842, 0.367220, 0.254952, 0.182860, 0.133060, 0.097244, 0.070868, 0.051199, 0.036466, 0.025435, 0.017157,
  exponents(:,2) =  0.437039, 0.303223, 0.217211, 0.158560, 0.116878, 0.086332, 0.063474, 0.046155, 0.032943, 0.022774,
  exponents(:,3) =  0.409883, 0.291455, 0.213760, 0.159587, 0.120199, 0.090678, 0.068084, 0.050566, 0.036876, 0.026055,
  exponents(:,4) =  0.340244, 0.243185, 0.179271, 0.134479, 0.101646, 0.076790, 0.057555, 0.042443, 0.030365,
  exponents(:,5) =  0.225709, 0.167874, 0.128505, 0.099316, 0.076719, 0.058780, 0.044303, 0.032392,
  exponents(:,6) =  0.118188, 0.085985, 0.062824, 0.045463, 0.032349, 0.022480, 0.015073,
  exponents(:,7) =  0.202350, 0.124177, 0.079606, 0.050353, 0.030895, 0.018109
/
&process_control
 max_ijrs_size = 35000.0
/
```

Namelist `pco_data`
-----------------------

* **max_PCO_l I 1 (-1)**

   The minimum partial wave for the PCOs which the PCO generation parameter arrays have as a lower index.

* **min_PCO_l I 1 (-1)**

   The maximum partial wave for the PCOs which the PCO exponent generation parameter arrays have as a upper index.

* **PCO_alpha0 R (min_PCO_l:max_PCO_l) (-1.0)**

   The alpha_0 PCO exponent generation parameter, as in formula alpha = alpha_0 * beta ** (i-1), where alpha is the exponent of the ith PCO shell, the value of alpha_0 should be < 1 and > 0.

* **PCO_beta R (min_PCO_l:max_PCO_l) (-1.0)**

   The beta PCO exponent generation parameter, as in formula alpha = alpha_0 * beta ** (i-1), where alpha is the exponent of the ith PCO shell, the value of beta should be > 1.

* **num_PCOs I (min_PCO_l:max_PCO_l) (-1)**

   The number of PCO shells/exponents per partial wave, gives i = 1, num_PCOs(l), as in formula alpha = alpha_0 * beta ** (i-1), where alpha is the exponent of the ith PCO shell.

* **PCO_gto_thrs R (min_PCO_l:max_PCO_l) (-1.0)**

   The threshold for eliminating continuum GTOs exponents too close to PCO exponents per partial wave, if not given default behaviour is to exclude any continuum GTO exponents closer to PCO exponents than alpha_0*(beta-1.0).

* **PCO_del_thrs R (8) (-1.0)**

   The threshold for deletion of PCO orbitals per irreduciable representation symmetry.

