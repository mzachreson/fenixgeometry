!***fgwrite.f90***



module write_mod

	implicit none
    contains


    !***********************
    ! Writing all the files:
    ! collstuff, parameters, spoly, masksource, maskcells, mcl, maskstuff

    subroutine write_geom_files
        use geombuilder_mod
        implicit none

        !Writing all the files:
        call write_collstuff
        call write_parameters
        call write_spoly
		call write_maskcells_and_mcl
        call write_maskstuff
		call write_samplingcells


    end subroutine write_geom_files








    !--write "collstuff.dat"
    subroutine write_collstuff
        use data_mod
        implicit none


        integer collstuff_unit
        integer i, j, k
        parameter( collstuff_unit = 201 )

        real(8) zero

        zero = 0.0

        write(*,*) "Writing 'collstuff.dat'..."
        open(unit = collstuff_unit, file="geom/collstuff.dat", status="replace")

        do i = 1,num_proc_r
            do j = 1,num_proc_z

                do k = 1,processors(i,j)%nc

                    write(collstuff_unit,"(5(1x,1pe16.8),i5,2(1x,1pe16.8),i5,2(1x,1pe16.8))")  processors(i,j)%collcells(k)%rmin * drfine, &
                                                                processors(i,j)%collcells(k)%rmax * drfine, &
                                                                processors(i,j)%collcells(k)%zmin * dzfine, &
                                                                processors(i,j)%collcells(k)%zmax * dzfine, &
                                                                processors(i,j)%collcells(k)%volume, &
                                                                processors(i,j)%collcells(k)%mask, &
                                                                zero, zero, 0, zero, zero     ! adding on a couple of zeros for Electric Field R, Electric Field Z
                                                                               ! These are modified later by "Ebuilder" if there is an added field.

                end do
            end do
        end do

        close(unit = collstuff_unit)

    end subroutine write_collstuff



    !--write "parameters.dat"
    subroutine write_parameters
        use data_mod
        implicit none

        integer i,j
        integer param_unit
        parameter(param_unit=202)

        integer :: nc_sum = 0

        write(*,*) "Writing 'parameters.dat'..."
        open(unit=param_unit, file='./geom/parameters.dat',status='replace')

        do i = 1,num_proc_r
            do j = 1,num_proc_z
                nc_sum = nc_sum + processors(i,j)%nc
            end do
        end do

        write(param_unit,*) nc_sum
        write(param_unit,*) (rmax-rmin)  !nrfine - number of drfine widths in rmax-rmin
        write(param_unit,*) (zmax-zmin)  !nzfine - number of dzfine widths in zmax-zmin
        write(param_unit,*) npoly
        write(param_unit,*) totnumpoints ! Number of mask codes
        write(param_unit,*) drfine
        write(param_unit,*) dzfine
        !write(param_unit,*) rmin * drfine  !we don't write this, in this version of FENIX (but I think we should in the future)
        write(param_unit,*) rmax * drfine
        write(param_unit,*) zmin * dzfine
        write(param_unit,*) zmax * dzfine
        

        close(unit=param_unit)
        

    end subroutine write_parameters




    !--write "spoly.dat"
    subroutine write_spoly
        use data_mod
        implicit none

        integer i, j
        integer spoly_unit
        parameter(spoly_unit=303)

        write(*,*) "Writing 'spoly.dat'..."
        open(unit=spoly_unit, file='./geom/spoly.dat',status='replace')

        do i = 1,npoly
            write(spoly_unit,*) polygons(i)%numpoints, polygons(i)%polytype, polygons(i)%reflect
            do j = 1,polygons(i)%numpoints
                write(spoly_unit,*) polygons(i)%points(j)%z*dzfine, polygons(i)%points(j)%r*drfine ! These are written FLIPPED..(Z,R) instead of (R,Z)
            end do
        end do

        close(unit=spoly_unit)

    end subroutine write_spoly




    !--write "maskcells.dat"
	! not actually read in by fenix.  I'm not sure what it's for, actually.
    subroutine write_maskcells_and_mcl
     use data_mod
        implicit none
		integer maskcells_unit, masknumber, i, j, k, line, linestart, lineend
		parameter(maskcells_unit = 404)

		integer mcl_unit
		parameter(mcl_unit = 505)

		line = 0

        write(*,*) "Writing 'maskcells.dat' and 'mcl.dat'..."
		open(unit=maskcells_unit, file='./geom/maskcells.dat',status='replace')
		open(unit=mcl_unit, file='./geom/mcl.dat',status='replace')

		! for each non-zero mask number, we loop over each
		! collision cell and print the indicies of the collision cells
		! with that mask number.
		do masknumber = 1,maxmaskcode
			linestart = line + 1
			do i = 1,num_proc_r
				do j = 1,num_proc_z

					do k = 1,processors(i,j)%nc
						if( processors(i,j)%collcells(k)%mask == masknumber ) then
							line = line + 1
							write(maskcells_unit,*) processors(i,j)%collcells(k)%psudoindex
						end if
					end do

				end do
			end do
			lineend = line

			! if there's at least one line of a particular 
			! mask code, then lineend >= linestart
			if( lineend >= linestart ) then
				write(mcl_unit,*) linestart, lineend
			else
				write(mcl_unit,*) 0, 0
			end if
		end do

		close(unit=maskcells_unit)
		close(unit=mcl_unit)

    end subroutine write_maskcells_and_mcl





    !--write "maskstuff.dat"
    subroutine write_maskstuff
        use data_mod
        implicit none

		integer maskstuff_unit, i
		parameter(maskstuff_unit = 606)

		open(unit=maskstuff_unit, file='./geom/maskstuff.dat',status='replace')
        write(*,*) "Writing 'maskstuff.dat'..."

		do i = 1,totnumpoints
			write(maskstuff_unit,"(1x,i5,2x,f6.4,11(1x,1pe16.8))")  masks(i)%code, masks(i)%reflect, &
										    masks(i)%s0z, masks(i)%s0r, &
		 								   masks(i)%w1z, masks(i)%w1r, &
										masks(i)%g1z, masks(i)%g1r, &
										masks(i)%w2z, masks(i)%w2r, &
										masks(i)%g2z, masks(i)%g2r, &
										masks(i)%temperature
		end do

		close(unit=maskstuff_unit)

    end subroutine write_maskstuff




    !**********
    ! A new one!  I call it "samplingcells.dat"
    ! with volume-corrected sampling cells.
	! format:   <i>   <j>   <volume>  <(bool) volume was corrected>
	subroutine write_samplingcells
		use data_mod
		use volume_calc_mod
		implicit none
		integer i,j
		real(8) volume, bottom, top, left, right
		logical vol_was_corrected
		integer nclr, nclz
		integer scl_unit, counter,flag
		parameter(scl_unit=707)

        write(*,*) "Writing 'samplingcells.dat'..."
		open(unit=scl_unit, file='./geom/samplingcells.dat',status='replace')

		nclr = (rmax-rmin) / s_cell_dr
		nclz = (zmax-zmin) / s_cell_dz

		counter = 0
		do i = 1,nclr
			do j = 1,nclz
				counter = counter + 1

				bottom = (i-1)*s_cell_dr
				top = (i)*s_cell_dr
				left = (j-1)*s_cell_dz
				right = (j)*s_cell_dz

                                flag=1 ! sampling cell volume correction
				volume = corrected_cell_volume( bottom, top, left, right, &
                                         drfine, dzfine, polygons, npoly, vol_was_corrected, counter, flag)

				write(scl_unit,*) i, j, "  ", volume, "  ", vol_was_corrected
			end do
		end do

		close(unit=scl_unit)

	end subroutine write_samplingcells


end module write_mod


