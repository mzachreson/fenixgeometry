!**** fgtypes.f90 ****





!--types_mod--
module types_mod
    implicit none
	
	! the findcell array is made up of these
	! so we know not only the index of the cell
	! but which processor it's in
	type index_type
		integer cl_index
		integer proc_i
		integer proc_j
	end type index_type


    !--a point, in drfine/dzfine units (even though it's a real(8))
    type point_type
        real(8) r
        real(8) z
    end type point_type


	! mask type -- this is what polygons get converted into.
	type mask_type
		integer code
		real(8) reflect

		real(8) s0r, s0z
		real(8) w1r, w1z, g1r, g1z
		real(8) w2r, w2z, g2r, g2z

		real(8) temperature
	end type mask_type


    !--a polygon
    type polygon_type
        integer numpoints
        integer polytype
        real(8) reflect
        real(8) temperature
        type(point_type), pointer, dimension(:) :: points
    end type polygon_type


    !--a region in Z in which cirtain
    !  characteristics are defined, which are
    !  unique to that region
    type zregion_type
        integer zmin, zmax
        real(8) smallest_cc

	    integer num_special_r
	    integer, allocatable, dimension(:) :: special_r
		integer, allocatable, dimension(:) :: spec_r_marker
    end type zregion_type


    !--a collision cell type
    type collcell_type
        integer rmin, rmax ! in drfine/dzfine units
        integer zmin, zmax ! in drfine/dzfine units
        real(8) volume ! in meters^3
        integer mask

		! this boolean specifies whether this cell's volume
		! has been corrected because of intersection with metal. (for debugging purposes)
		logical volcorrected

		! "psudoindex" is the line number you will find this cell in the "collcells.dat" file
		integer psudoindex

    end type collcell_type


    !--processor_type
    type processor_type
        integer glob_rmin, glob_rmax ! in drfine/dzfine units
        integer glob_zmin, glob_zmax ! in drfine/dzfine units

        integer rmin, rmax ! in drfine/dzfine units
        integer zmin, zmax ! in drfine/dzfine units
        real(8) drfine, dzfine

        integer nreg
        type(zregion_type), pointer, dimension(:) :: zregions

        integer num_special_r
        integer, pointer, dimension(:) :: special_r

        integer nc  ! number of collision cells
		integer ccarraysize ! size of the collcells array
        type(collcell_type), pointer, dimension(:) :: collcells

        ! every processor_type has a copy of the entire
        ! simulation's polygons.  This is used
        ! to calculate the mask numbers, and corrected volumes.
        integer npoly
        type(polygon_type), pointer, dimension(:) :: polygons

    end type processor_type


end module types_mod




