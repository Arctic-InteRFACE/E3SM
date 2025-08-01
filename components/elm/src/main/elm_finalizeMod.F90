module elm_finalizeMod

  !-----------------------------------------------------------------------
  ! Performs land model cleanup
  !
  !
  use SoilLittVertTranspMod, only : cleanupLitterTransportList

  implicit none
  save
  public   ! By default everything is public
  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------
  !
  public :: final
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine final( )
    !
    ! !DESCRIPTION:
    ! Finalize land surface model
    !
#ifdef USE_PETSC_LIB
#include <petsc/finclude/petsc.h>
#endif
    ! !USES:
    use elm_varctl             , only : use_vsfm, use_cn
#ifdef USE_PETSC_LIB
    use petscsys
#endif
    !
    ! !ARGUMENTS
    implicit none
    !

#ifdef USE_PETSC_LIB
    PetscErrorCode        :: ierr

    if (use_vsfm) then
       call PetscFinalize(ierr)
    endif
#endif
    if (use_cn) then
       call cleanupLitterTransportList()
    endif

  end subroutine final

end module elm_finalizeMod
