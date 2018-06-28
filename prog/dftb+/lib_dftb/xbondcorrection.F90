!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Halogen bond correction for DFTB3 - additional repulsive potential applied to Cl,Br,I - O,N
!> pairs. See http://dx.doi.org/10.1021/ct5009137 for details
module xbondcorrection
  use accuracy
  use constants
  use vdwdata
  use io
  use message, only : warning
  implicit none
  private

  public :: XBCorr, XBCorr_init

  !> Internal data of the X-bond correction
  type :: XBCorr
    private

    ! Matrix of pairwise parameters between species
    real(dp), allocatable :: xbPairParams(:,:)

    ! Matrix of sums of covalent radii
    real(dp), allocatable :: xbSpeciesRadii(:)

    !> Calculated correction to energy
    real(dp) :: xbEnergy

  contains

    procedure :: updateCoords
    procedure :: calcXBEnergy

  end type XBCorr

  !> Global parameters (original values from http://dx.doi.org/10.1021/ct5009137
  !> in Angstrom and kcal/mol unit system converted to a.u.)
  real(dp), parameter :: xbGlobalC1 = 7.761_dp * kcal_mol__Hartree
  real(dp), parameter :: xbGlobalC3 = 4.518_dp ! Dimensionless
  real(dp), parameter :: xbGlobalC2 = 0.050_dp * (1.0_dp / AA__Bohr)**xbGlobalC3
  real(dp), parameter :: xbRCut0 = 0.7_dp ! Dimensionless
  real(dp), parameter :: xbRCut1 = 0.8_dp ! Dimensionless
 
  !> Pairwise parameters (in Angstrom, converted to Bohrs)
  real(dp), parameter :: xbParamOCl = 1.237_dp * AA__Bohr
  real(dp), parameter :: xbParamOBr = 1.099_dp * AA__Bohr
  real(dp), parameter :: xbParamOI = 1.313_dp * AA__Bohr
  real(dp), parameter :: xbParamNCl = 1.526_dp * AA__Bohr
  real(dp), parameter :: xbParamNBr = 1.349_dp * AA__Bohr
  real(dp), parameter :: xbParamNI = 1.521_dp * AA__Bohr

  !> vdW radii used (in Angstrom, converted to Bohrs)
  real(dp), parameter :: xbRadiusO = 1.52_dp * AA__Bohr
  real(dp), parameter :: xbRadiusN = 1.55_dp * AA__Bohr
  real(dp), parameter :: xbRadiusCl = 1.75_dp * AA__Bohr
  real(dp), parameter :: xbRadiusBr = 1.85_dp * AA__Bohr
  real(dp), parameter :: xbRadiusI = 1.98_dp * AA__Bohr


contains

  subroutine XBCorr_init(this, speciesNames)

    !> Initialised instance at return.
    type(XBCorr), intent(out) :: this

    !> Names of the species
    character(mc), allocatable, intent(in) :: speciesNames(:)

    integer :: nSpecies
    integer :: iSp1, iSp2
    character(mc) :: spName1, spName2

    ! Allocate matrix of parameters
    nSpecies = size(speciesNames)
    allocate(this%xbPairParams(nSpecies, nSpecies))
    allocate(this%xbSpeciesRadii(nSpecies))
    ! Fill with -1 (no correction applied)
    this%xbPairParams(:,:) = -1.0_dp
    ! Copy pairwise parameters
    do iSp1 = 1, size(speciesNames)
      do iSp2 = 1, size(speciesNames)
        spName1 = speciesNames(iSp1)
        spName2 = speciesNames(iSp2)

        if (count([spName1, spName2] == "O") == 1) then
          ! Pairs with O
          if (count([spName1, spName2] == "Cl") == 1) then
            this%xbPairParams(iSp1,iSp2) = xbParamOCl
            this%xbPairParams(iSp2,iSp1) = xbParamOCl
          else if (count([spName1, spName2] == "Br") == 1) then
            this%xbPairParams(iSp1,iSp2) = xbParamOBr
            this%xbPairParams(iSp2,iSp1) = xbParamOBr
          else if (count([spName1, spName2] == "I") == 1) then
            this%xbPairParams(iSp1,iSp2) = xbParamOI
            this%xbPairParams(iSp2,iSp1) = xbParamOI
          end if
        else if (count([spName1, spName2] == "N") == 1) then
          ! Pairs with N
          if (count([spName1, spName2] == "Cl") == 1) then
            this%xbPairParams(iSp1,iSp2) = xbParamNCl
            this%xbPairParams(iSp2,iSp1) = xbParamNCl
          else if (count([spName1, spName2] == "Br") == 1) then
            this%xbPairParams(iSp1,iSp2) = xbParamNBr
            this%xbPairParams(iSp2,iSp1) = xbParamNBr
          else if (count([spName1, spName2] == "I") == 1) then
            this%xbPairParams(iSp1,iSp2) = xbParamNI
            this%xbPairParams(iSp2,iSp1) = xbParamNI
          end if
        end if
      end do
    end do
    ! Copy species radii
    write (stdout, *) "Radii:"
    do iSp1 = 1, size(speciesNames)
      write (stdout, *) speciesNames(iSp1)
      call getVdwData(speciesNames(iSp1), this%xbSpeciesRadii(iSp1))
      write (stdout, *) this%xbSpeciesRadii(iSp1)

      select case (speciesNames(iSp1))
      case ("O")
        this%xbSpeciesRadii(iSp1) = xbRadiusO
      case ("N")
        this%xbSpeciesRadii(iSp1) = xbRadiusN
      case ("Cl")
        this%xbSpeciesRadii(iSp1) = xbRadiusCl
      case ("Br")
        this%xbSpeciesRadii(iSp1) = xbRadiusBr
      case ("I")
        this%xbSpeciesRadii(iSp1) = xbRadiusI
      end select
      write (stdout, *) this%xbSpeciesRadii(iSp1)
    end do

    ! Printing of parameters, for debugging only
    if (.false.) then
      write (stdout, *) "Initializing X-bond correction, parameter matrix:"
      do iSp1 = 1, size(speciesNames)
        do iSp2 = 1, size(speciesNames)
          spName1 = speciesNames(iSp1)
          spName2 = speciesNames(iSp2)
          write (stdout, *) trim(spName1) // "-" // trim(spName2) // " = ", this%xbPairParams(iSp1,iSp2)
        end do
      end do
      write (stdout, *)
    end if

  end subroutine XBCorr_init


  !> Upon change of coordinates, calculate the correction
  subroutine updateCoords(this, coords, species)
    !> instance of the correction
    class(XBCorr), intent(inout) :: this

    !> Current coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Species of atoms
    integer, allocatable, intent(in) :: species(:)

    integer :: nAtom, i, j
    real(dp) :: vect(3), r, r_vdw, energy, r0, r1, sw, x, erep, erep0

    ! Get no. of atoms
    nAtom = size(coords, dim=2)
   
    ! Iterate all pairs of atoms 
    energy = 0.0_dp
    do i = 1, nAtom
      do j = 1, i-1
        ! Skip pairs to which correction is not applied
        if (this%xbPairParams(species(i),species(j)) == -1.0_dp) then
          cycle
        end if
        ! Calculate interatomic vector and distance
        vect(:) = coords(:,i) - coords(:,j)
        r = sqrt(sum(vect**2))
        ! Sum of vdW radii
        r_vdw = this%xbSpeciesRadii(species(i)) + this%xbSpeciesRadii(species(j))
        ! Switching function
        r0 = xbRCut0 * r_vdw
        r1 = xbRCut1 * r_vdw
        if (r <= r0) then
          sw = 0.0_dp
        else if (r >= r1) then
          sw = 1.0_dp
        else
          x = (r - r0) / (r1 - r0)
          sw = -20.0_dp * x**7 + 70.0_dp * x**6 - 84.0_dp * x**5 + 35.0_dp * x**4
        end if
        ! Calculate current and cutoff values
        ! the max(...,0) is added to the original implementation to allow secure calculation at short distances
        erep = xbGlobalC1 * exp(-xbGlobalC2 * max(r - this%xbPairParams(species(i),species(j)),0.0_dp)**xbGlobalC3)
        erep0 = xbGlobalC1 * exp(-xbGlobalC2 * max(r0 - this%xbPairParams(species(i),species(j)),0.0_dp)**xbGlobalC3)
        ! Apply the switching function
        erep = erep * sw + erep0 * (1.0_dp - sw)
        energy = energy + erep
      end do
    end do
    ! Save calculated energy
    this%xbEnergy = energy
  end subroutine updateCoords

  !> Get the calculated energy and apply it to the results
  subroutine calcXBEnergy(this, energy)
    !> instance of the correction
    class(XBCorr), intent(inout) :: this
    !> Energy to be corrected
    real(dp), intent(inout) :: energy

    energy = energy + this%xbEnergy
  end subroutine calcXBEnergy
end module xbondcorrection
