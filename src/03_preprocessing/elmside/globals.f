!***********************************************************
!         A module to contain the variables
!         so that they can be used globally    
!         This saves on sub call arguments.
!***********************************************************

module globals
    implicit none

    integer, parameter  :: maxneighbours2D = 20
    integer, parameter  :: maxneighbours3D = 150

! to store file data
    integer                            :: ndimensions
    integer                            :: nvertices   ! number of match points
    integer, allocatable, dimension(:) :: vertices    ! list of points to match to

    integer                            :: nelems
    integer, allocatable, dimension(:,:) :: v ! vertex ID per elem.
    integer                            :: d1,d2,d3    ! dummies


    integer                            :: nallvertices ! number of total points in mesh
    double precision, allocatable, dimension(:,:) :: coords
    double precision                   :: xmin, xmax ! domain extrema
    double precision                   :: ymin, ymax ! domain extrema
    double precision                   :: zmin, zmax ! domain extrema

    integer                            :: numloop
    double precision, allocatable, dimension(:,:) :: loopcoords

    integer, allocatable, dimension(:) :: nneighbours
    integer, allocatable, dimension(:,:) :: neighbourIDs

    integer, dimension(10)              :: neighbours

! orthogonal vector direction
    double precision                   :: ax, ay, az          ! anchor point
    double precision                   :: ux, uy, uz          ! vector to span up the base
    double precision                   :: vx, vy, vz          ! vector to span up the base, only for 3D
    double precision                   :: ox, oy, oz, normo   ! orthogonal vector to the base, and its norm
    integer                            :: pm_sign             ! sign of the ortho vec +/- 1
    double precision                   :: eps

! flags from command line options
    logical                            :: printcoords
    logical                            :: printsides
! attachment mode. Give elements that are attached with only their side (sideonly == .true.),
! or even with a point (sideonly == .false.)
    logical                            :: sideonly
    logical                            :: opn
    logical                            :: printweights
    logical                            :: invertsign
    logical                            :: endpoints

! filenames passed in arguments
    character (len=128)                :: tecindatpartfnps  ! optional
    character (len=128)                :: tecindatpartfelm  ! mandatory
    character (len=128)                :: vertexlist        ! mandatory
    character (len=128)                :: loopfile          ! optional

! level of output
    integer                            :: iecho

! algorithm used
    integer                            :: algorithm

    integer, parameter                 :: stringLength = 500

    character(len=stringLength)        :: output(4)



! subs

    public :: freeglobals

    contains

    subroutine freeglobals

        implicit none

        if (allocated(loopcoords)) then
            deallocate(loopcoords)
        endif
        deallocate(coords)
        deallocate(vertices)
        deallocate(v)
        deallocate(neighbourIDs)
        deallocate(nneighbours)
    end subroutine

end module

