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

    ! Van der Waals radii of the species
    real(dp), allocatable :: xbSpeciesRadii(:)

    !> Calculated correction to energy
    real(dp) :: xbEnergy

    !> Derivatives of the correction energy
    real(dp), allocatable :: xbDerivs(:,:)

  contains

    procedure :: updateCoords
    procedure :: calcXBEnergy
    procedure :: addGradients

  end type XBCorr

  !> Parameters (original values from http://dx.doi.org/10.1021/ct5009137
  !> in Angstrom and kcal/mol unit system converted to a.u.)

  !> Global parameters
  real(dp), parameter :: xbGlobalC1 = 7.761_dp * kcal_mol__Hartree
  real(dp), parameter :: xbGlobalC3 = 4.518_dp ! Dimensionless
  real(dp), parameter :: xbGlobalC2 = 0.050_dp * (1.0_dp / AA__Bohr)**xbGlobalC3
  real(dp), parameter :: xbRCut0 = 0.7_dp ! Dimensionless
  real(dp), parameter :: xbRCut1 = 0.8_dp ! Dimensionless
 
  !> Pairwise parameters
  real(dp), parameter :: xbParamOCl = 1.237_dp * AA__Bohr
  real(dp), parameter :: xbParamOBr = 1.099_dp * AA__Bohr
  real(dp), parameter :: xbParamOI = 1.313_dp * AA__Bohr
  real(dp), parameter :: xbParamNCl = 1.526_dp * AA__Bohr
  real(dp), parameter :: xbParamNBr = 1.349_dp * AA__Bohr
  real(dp), parameter :: xbParamNI = 1.521_dp * AA__Bohr

contains

  subroutine XBCorr_init(this, speciesNames)

    !> Initialised instance at return.
    type(XBCorr), intent(out) :: this

    !> Names of the species
    character(mc), allocatable, intent(in) :: speciesNames(:)

    integer :: nSpecies
    integer :: iSp1, iSp2
    character(mc) :: spName1, spName2

    ! Build a matrix of pairwise parameters
    nSpecies = size(speciesNames)
    allocate(this%xbPairParams(nSpecies, nSpecies))
    allocate(this%xbSpeciesRadii(nSpecies))
    ! Fill with -1 (no correction applied)
    this%xbPairParams(:,:) = -1.0_dp
    ! Copy pairwise parameters for current system
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
    ! I have checked that the radii in DFTB+'s database match my radii
    ! for all the elements used
    do iSp1 = 1, size(speciesNames)
      call getVdwData(speciesNames(iSp1), this%xbSpeciesRadii(iSp1))
    end do

  end subroutine XBCorr_init


  !> Upon change of coordinates, calculate the correction
  subroutine updateCoords(this, coords, species)
    !> instance of the correction
    class(XBCorr), intent(inout) :: this

    !> Current coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Species of atoms
    integer, allocatable, intent(in) :: species(:)

    integer :: nAtom, i, j, k
    real(dp) :: vect(3), r, r_vdw, energy, r0, r1, sw, x, d, erep, erep0
    real(dp) :: dsw_dr, de_dr, cart_deriv(3)

    ! Get no. of atoms
    nAtom = size(coords, dim=2)

    ! Alllocate and reset the derivatives
    if (.not. allocated(this%xbDerivs)) then
      allocate(this%xbDerivs(3,nAtom))
    end if
    this%xbDerivs(:,:) = 0.0_dp
   
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
        d = this%xbPairParams(species(i),species(j))
        erep = xbGlobalC1 * exp(-xbGlobalC2 * (r - d)**xbGlobalC3)
        erep0 = xbGlobalC1 * exp(-xbGlobalC2 * (r0 - d)**xbGlobalC3)
        ! Apply the switching function
        erep = erep * sw + erep0 * (1.0_dp - sw)
        energy = energy + erep

        ! Gradient
        dsw_dr = 0.0_dp
        if (r > r0 .and. r < r1) then
          dsw_dr = -140.0_dp * x**6 + 420.0_dp * x**5 - 420.0_dp * x**4 + 140.0_dp * x**3
        end if
        de_dr = -1.0_dp * xbGlobalC2 * xbGlobalC3 * (r - d)**(xbGlobalC3 - 1.0_dp) * erep
        de_dr = sw * de_dr + dsw_dr * erep + dsw_dr * erep0
        ! Convert gradient co cartesian and add it to the complete gradient
        cart_deriv = vect * (de_dr/r)
        do k = 1, 3
          this%xbDerivs(k,i) = this%xbDerivs(k,i) + cart_deriv(k)
          this%xbDerivs(k,j) = this%xbDerivs(k,j) - cart_deriv(k)
        end do
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

  subroutine addGradients(this, derivs)
    !> instance of the correction
    class(XBCorr), intent(inout) :: this

    !> Derivatives to be modified
    real(dp), intent(inout) :: derivs(:,:)

    derivs = derivs + this%xbDerivs

  end subroutine addGradients

end module xbondcorrection
