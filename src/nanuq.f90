PROGRAM nanuq
   !!======================================================================
   !!                     ***  PROGRAM nanuq  ***
   !!
   !! ** Purpose :   encapsulate nanuq_gcm so that it can also be called
   !!              together with the linear tangent and adjoint models
   !!======================================================================
   !! History :   OPA  ! 2001-02  (M. Imbard, A. Weaver)  Original code
   !!   SI3      1.0  ! 2003-10  (G. Madec) F90
   !!----------------------------------------------------------------------
   USE nanuqgcm   ! NANUQ system   (nanuq_gcm routine)
   !!----------------------------------------------------------------------
   !! NANUQ 0.1 beta, Brodeau (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
   !
   CALL nanuq_gcm           ! NANUQ direct code
   !
   !!======================================================================
END PROGRAM nanuq
