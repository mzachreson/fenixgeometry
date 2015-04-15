!**** fgdata.f90 ******



!--data_mod--
! Here all the data obtained from "input_data.txt" is stored
module data_mod
    use types_mod
    implicit none
    save

    ! data fread in from the file:
    real(8) Nef, vth
    
    integer nsteps
    real(8) tau

    logical use_restart, single_proc_mode
    integer output_skip, restart_skip
    real(8) drfine, dzfine
    integer s_cell_dr, s_cell_dz
    integer rmin, rmax, zmin, zmax
    integer num_proc_r, num_proc_z
    integer, allocatable, dimension(:) :: comp_reg_r, comp_reg_z

    integer nreg
    type(zregion_type), allocatable, dimension(:) :: zregions

!These are being moved into the "zregions" data type.
!    integer num_special_r
!    integer, allocatable, dimension(:) :: special_r

    integer npoly
    type(polygon_type), allocatable, dimension(:) :: polygons


    ! used in geombuilder_mod
    type(processor_type), allocatable, dimension(:,:) :: processors
	type(index_type), allocatable, dimension(:,:) :: findcell

	integer maxmaskcode
	integer totnumpoints ! this is the total number of geometry points
						 ! in the geometry, and also is the exact number of masks created.
	type(mask_type), allocatable, dimension(:) :: masks







    contains

    subroutine read_input_data()
        use types_mod
        implicit none

        !note: everything INTERNALLY is in INTEGER UNITS OF dzfine, drfine
        ! except: Nef, vth
        !         dzfine, drfine, smallest_cc

        ! input file unit number:
        integer input_unit
        parameter(input_unit=543)

        ! used for do loops:
        integer i, j

        ! used because some things are entered as integers, in units
        ! of either dzfine/drfine or s_cell_dr/s_cell_dz
        integer sampling_int_R, sampling_int_Z
        integer rmin_int, rmax_int, zmin_int, zmax_int
        integer, allocatable, dimension(:) :: zregright_int
        integer pointr, pointz


        open(unit=input_unit, file="input_data.txt", status="old")

        read(input_unit,*) Nef, vth
        read(input_unit,*) nsteps, tau
        read(input_unit,*) use_restart, single_proc_mode
        read(input_unit,*) output_skip, restart_skip
        !Skip two lines of input_data.txt.  These lines have information 
        !That is only needed by the simulation
        read(input_unit,*) 
        read(input_unit,*)
        read(input_unit,*) drfine, dzfine
        read(input_unit,*) s_cell_dr, s_cell_dz

        read(input_unit,*) rmin, rmax, zmin, zmax
        ! convert to units of dzfine, drfine
        rmin = rmin * s_cell_dr
        rmax = rmax * s_cell_dr
        zmin = zmin * s_cell_dz
        zmax = zmax * s_cell_dz

        !************ Parallel stuff:

        read(input_unit,*) num_proc_r, num_proc_z

        allocate(comp_reg_r(num_proc_r+1))
        allocate(comp_reg_z(num_proc_z+1))

        comp_reg_r(1) = 0
        comp_reg_z(1) = 0

        do i = 2,num_proc_r+1
            read(input_unit,*) comp_reg_r(i)
            ! converting to units of drfine
            comp_reg_r(i) = comp_reg_r(i) * s_cell_dr
			!write(*,*) comp_reg_r(i)
        end do

        do i = 2,num_proc_z+1
            read(input_unit,*) comp_reg_z(i)
            ! converting to units of dzfine
            comp_reg_z(i) = comp_reg_z(i) * s_cell_dz
			!write(*,*) comp_reg_z(i)
        end do
        !************ zregion stuff:

        read(input_unit,*) nreg
        allocate(zregions(nreg))
        allocate(zregright_int(nreg))

        ! zrightreg:
        read(input_unit,*) (zregright_int(i),i=1,nreg)
        zregions(1)%zmin = zmin
        do i = 1,nreg
            if( i > 1 ) zregions(i)%zmin = zregions(i-1)%zmax 
            zregions(i)%zmax = zregright_int(i) * s_cell_dz
        end do

        ! smallest_cc:
        read(input_unit,*) (zregions(i)%smallest_cc,i=1,nreg)

        ! special_r:
        read(input_unit,*) (zregions(i)%num_special_r, i=1,nreg)
		do i = 1,nreg
        	if( zregions(i)%num_special_r > 0 ) then
            	allocate(zregions(i)%special_r(zregions(i)%num_special_r))
            	read(input_unit,*) (zregions(i)%special_r(j), j=1,zregions(i)%num_special_r)
				do j = 1, zregions(i)%num_special_r
					! converting to drfine units
					zregions(i)%special_r(j) = zregions(i)%special_r(j) * s_cell_dr
				end do
        	end if
		end do
        !*********** Polygon stuff:

        read(input_unit,*) npoly
        allocate(polygons(npoly))
        do i = 1,npoly
            read(input_unit,*) polygons(i)%numpoints, polygons(i)%polytype, polygons(i)%reflect, polygons(i)%temperature
            allocate(polygons(i)%points(polygons(i)%numpoints))
            do j = 1,polygons(i)%numpoints
                read(input_unit,*) pointr, pointz

                ! rounding point r to the nearest integer of drfine
                if( (pointr * s_cell_dr) - dble(int( pointr * s_cell_dr )) >= 0.5 ) then
                    polygons(i)%points(j)%r = int(pointr * s_cell_dr) + 1
                else
                    polygons(i)%points(j)%r = int(pointr * s_cell_dr)
                end if

                ! rounding point z to the nearest integer of dzfine
                if( (pointz * s_cell_dz) - dble(int( pointz * s_cell_dz)) >= 0.5 ) then
                    polygons(i)%points(j)%z = int(pointz * s_cell_dz) + 1
                else
                    polygons(i)%points(j)%z = int(pointz * s_cell_dz)
                end if
                
            end do
        end do

        close(unit=input_unit)


    end subroutine read_input_data





end module data_mod

