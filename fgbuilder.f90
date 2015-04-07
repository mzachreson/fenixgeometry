!***fgbuilder.f90***



!--geombuilder_mod--
! uses the data in "data_mod" to process
! the geometry into the proper form, so that
! everything is consistant.
module geombuilder_mod
	use volume_calc_mod
    use types_mod
    use data_mod
    implicit none

    contains

    !******** process_geom_data ****************
    !This subroutine has several purposes:
    !  1. Process all the "data_mod" data:
    !       Correct any part of it that's not quite consistant
    !       (such as shifting the geometry to fit the sampling cells or collision cells)
    !       Correct the volumes of collision cells and sampling cells which are partly inside metals.
    !  2. Load all the data into the same types of data structures as are used in FENIX, after reading in the geom files
    !  3. Write the data to the geometry files
    subroutine process_geom_data
		use mask_calc_mod
        use data_mod
        implicit none
        integer :: num_boundaries_inside = 0
        integer regsave(nreg)
        integer i,j,k,s,n, tot_nc, temp_int
		integer ipos, jpos
		type(collcell_type) tmpcl
		integer ijoint, iwall, imask
		real(8) paintstep

		
    
        ! creating the processors
        allocate( processors(num_proc_r, num_proc_z) )

        ! and organizing everything inside each of them
        tot_nc = 0
        do i = 1,num_proc_r
            do j = 1,num_proc_z

                write(*,*) ""
                write(*,*) "WORKING ON PROCESSOR ", i, j
                ! dimensions of the entire simulation region, for reference
                processors(i,j)%glob_rmin = rmin
                processors(i,j)%glob_rmax = rmax
                processors(i,j)%glob_zmin = zmin
                processors(i,j)%glob_zmax = zmax

                ! this processor's local dimensions
                processors(i,j)%rmin = comp_reg_r(i)
                processors(i,j)%rmax = comp_reg_r(i+1)
                processors(i,j)%zmin = comp_reg_z(j)
                processors(i,j)%zmax = comp_reg_z(j+1)

                ! for local reference, we give it a copy of drfine and dzfine
                processors(i,j)%drfine = drfine
                processors(i,j)%dzfine = dzfine

                ! Giving this processor a copy of the polygons of the geometry
                ! so that it can be used to calculate the mask numbers later.
                processors(i,j)%npoly = npoly
                allocate( processors(i,j)%polygons( npoly ) )
                do k = 1,npoly
                    processors(i,j)%polygons(k) = polygons(k)
                end do


                ! zregions -- if a zregion changes in the middle of
                !                a processor, then the processor will
                !                have more than one zregion

                ! first we run the loop once, to find out how big the array needs to be for this processor
                num_boundaries_inside = 0
                do k = 1, nreg

                    if(     ( zregions(k)%zmin <= processors(i,j)%zmin .and. zregions(k)%zmax <= processors(i,j)%zmin ) &
                       .or. ( zregions(k)%zmin >= processors(i,j)%zmax .and. zregions(k)%zmax >= processors(i,j)%zmax ) ) then
                        ! do nothing
                    else    
                        num_boundaries_inside = num_boundaries_inside + 1
						! the "regsave" array links which index in the whole region refers to which index in just the specific processor.
                        regsave(num_boundaries_inside) = k
                    end if
                end do

                ! check
                if( num_boundaries_inside < 1 ) then
                    write(*,*) "processor should have at least one region!!"
                    write(*,*) "n_b_inside:", num_boundaries_inside
                    write(*,*) "processor ", i, j
                    write(*,*) "proc zmin: ", processors(i,j)%zmin
                    write(*,*) "proc zmax: ", processors(i,j)%zmax
                    stop 667
                end if

                ! now we know how many regions are in the processor
                processors(i,j)%nreg = num_boundaries_inside
                allocate( processors(i,j)%zregions(num_boundaries_inside) )

                ! defining the zregion boundaries in the processor
                do n = 1,num_boundaries_inside
                    ! bottom boundary of the region
                    if( processors(i,j)%zmin >= zregions(regsave(n))%zmin ) then
                        processors(i,j)%zregions(n)%zmin = processors(i,j)%zmin
                    else
                        processors(i,j)%zregions(n)%zmin = zregions(regsave(n))%zmin
                    end if
                    
                    ! top boundary of the region
                    if( processors(i,j)%zmax <= zregions(regsave(n))%zmax ) then
                        processors(i,j)%zregions(n)%zmax = processors(i,j)%zmax
                    else
                        processors(i,j)%zregions(n)%zmax = zregions(regsave(n))%zmax
                    end if

                    !other information for the region
                    processors(i,j)%zregions(n)%smallest_cc = zregions(regsave(n))%smallest_cc

					!special_r values
                	if( zregions(regsave(n))%num_special_r > 0 ) then
                    	! finding out how many special_rs are in this processor region
                    	do k = 1,zregions(regsave(n))%num_special_r
                    	    if( zregions(regsave(n))%special_r(k) > processors(i,j)%rmin &
											.and. zregions(regsave(n))%special_r(k) < processors(i,j)%rmax ) then

                    	        processors(i,j)%zregions(n)%num_special_r = processors(i,j)%zregions(n)%num_special_r + 1

                    	    end if
                    	end do
	
                    	! allocating the array
                    	allocate( processors(i,j)%zregions(n)%special_r( processors(i,j)%zregions(n)%num_special_r ) )


                    	! putting the special_r values in
                    	temp_int = 0
                    	do k = 1,zregions(regsave(n))%num_special_r
                    	    if( zregions(regsave(n))%special_r(k) > processors(i,j)%rmin &
											.and. zregions(regsave(n))%special_r(k) < processors(i,j)%rmax ) then

                    	        temp_int = temp_int + 1
                    	        processors(i,j)%zregions(n)%special_r(temp_int) = zregions(regsave(n))%special_r(k)

                    	    end if
                    	end do
                	end if

                end do


                ! now the processor is initialized enough to be sent to have
                ! it's collision cells made
				
                call make_coll_cells( processors(i,j), tot_nc )
                tot_nc = tot_nc + processors(i,j)%nc

            end do
        end do

        write(*,*) ""
        write(*,*) "Total Number of Collision Cells:", tot_nc
        write(*,*) ""





		!now we will moke a findcells array
		! findcell(i,j) is the index of the cell where (drfine,dzfine) of (i,j) is located.
		write(*,*) "Creating findcell array..."

		allocate( findcell( rmax-rmin, zmax-zmin ) )


		do i = 1,num_proc_r
            do j = 1,num_proc_z
				do k = 1,processors(i,j)%nc


					! here we assign values in the findcells array
					tmpcl = processors(i,j)%collcells(k)
					do ipos = tmpcl%rmin+1,tmpcl%rmax
						do jpos = tmpcl%zmin+1,tmpcl%zmax
							findcell(ipos,jpos)%cl_index = k
							findcell(ipos,jpos)%proc_i = i
							findcell(ipos,jpos)%proc_j = j
						end do
					end do


					!********
				end do
			end do
		end do
		write(*,*) "  ... Done creating findcell array."

		


		! This is just a check, to make sure all of the indicies of the find cell
		! array have information in them.
		write(*,*) "Checking the findcell array..."
		do i = 1,rmax-rmin
			do j = 1,zmax-zmin
				if( 	 findcell(i,j)%proc_i < 1 .or. findcell(i,j)%proc_i > num_proc_r &
					.or. findcell(i,j)%proc_j < 1 .or. findcell(i,j)%proc_j > num_proc_z &
					.or. findcell(i,j)%cl_index < 1 ) then

					write(*,*) ""
					write(*,*) "ERROR!!  Find Cell Array Has problems!!!!"
					write(*,*) "Here's an example:"
					write(*,*) num_proc_r, "num_proc_r"
					write(*,*) num_proc_z, "num_proc_z"
					write(*,*) rmax-rmin, "rmax-rmin"
					write(*,*) zmax-zmin, "zmax-zmin"
					write(*,*) i,j, "i,j"
					write(*,*) findcell(i,j)%proc_i, "findcell(i,j)%proc_i"
					write(*,*) findcell(i,j)%proc_j, "findcell(i,j)%proc_j"
					write(*,*) findcell(i,j)%cl_index, "findcell(i,j)%cl_index"

					write(*,*) "        Writing file (called 'badcells.txt') of coordinates of offending indicies..."
					open( unit = 314159, file='badcells.txt', status='new' )
					do k = 1,rmax-rmin
						do s = 1,zmax-zmin
							if( 	 findcell(k,s)%proc_i < 1 .or. findcell(k,s)%proc_i > num_proc_r &
								.or. findcell(k,s)%proc_j < 1 .or. findcell(k,s)%proc_j > num_proc_z &
								.or. findcell(k,s)%cl_index < 1 ) then

								write(314159,*) k, s
							end if
						end do
					end do
					close( unit = 314159 )
					stop 2345
					
				end if
			end do
		end do
		write(*,*) "  ... Done checking the findcell array."



		! assigning the mask numbers to the cells... (mask painting)


		totnumpoints = 0
		do i = 1,npoly
			totnumpoints = totnumpoints + polygons(i)%numpoints
		end do

		! allocate two times as many mask slots as polygons, in case
		! they're all odd or all even, or something.
		allocate( masks(totnumpoints) ) ! there's the same number of masks as polygons
		maxmaskcode = 0 	! initialized to zero.

		!initialize everything to zero
		do i = 1,totnumpoints
			masks(i)%code = 0
			masks(i)%reflect = 0.0
			masks(i)%s0r = 0.0
			masks(i)%s0z = 0.0
			masks(i)%w1r = 0.0
			masks(i)%w1z = 0.0
			masks(i)%g1r = 0.0
			masks(i)%g1z = 0.0
			masks(i)%w2r = 0.0
			masks(i)%w2z = 0.0
			masks(i)%g2r = 0.0
			masks(i)%g2z = 0.0
			masks(i)%temperature = 0.0
		end do

		iwall = 0
		ijoint = -1
		imask = 0
		paintstep = 8.0*vth*tau


		write(*,*) ""
		write(*,*) "Painting masks..."

		do i = 1,npoly
			! here we adapt my code to work with Dr. Spencer's subroutine "maskpaint"
			call maskpaint( polygons(i), paintstep, iwall, ijoint, imask )
		end do
		write(*,*) "   ... Done painting masks."


                !Now call a subroutine that assigns a negative integer mask code
                !to  every cell inside the polygons (inside metal) with a negative number.  
                !A mask of -1 corresponds to polygon #1, -2 to polygon #2, etc.
                call mask_metal(polygons)  !written in fgmaskcalc.f90



!		open( unit=109, file = "test_vol_output.txt", status="replace")
!		do i = 1,num_proc_r
!			do j = 1,num_proc_z
!				!write(109,*) processors(i,j)%rmin, processors(i,j)%rmax, processors(i,j)%zmin, processors(i,j)%zmax
!				do k = 1,processors(i,j)%nc
!					if( processors(i,j)%collcells(k)%volcorrected ) then
!						write(109,*) processors(i,j)%collcells(k)%rmin, processors(i,j)%collcells(k)%zmin
!					end if
!				end do
!			end do
!		end do
!		close(unit=109)

	end subroutine process_geom_data
















    ! --- make_coll_cells ---- makes the collision cells
    ! And mask numbers, too.
    ! I'm going to make this particular subroutine "modular"
    subroutine make_coll_cells( node, cellcounter_prev_sum )
        use types_mod
		use volume_calc_mod
        implicit none
        type(processor_type) node
		integer cellcounter_prev_sum  ! the sum of all the cells of processors gone before


        !*****
        integer :: rmin
        integer :: rmax

        real(8) :: drfine
        real(8) :: dzfine
        integer :: nreg

        real(8) :: smallest_cc(node%nreg)
        integer :: regzmin(node%nreg)
        integer :: regzmax(node%nreg)

        integer :: reg_width(node%nreg)
        integer :: reg_height !the same for all regions
        integer :: numwide(node%nreg)
        integer :: cell_width(node%nreg)
		real(8) :: cell_width_real(node%nreg)
		real(8) :: width_remnant
		integer :: rightedge, leftedge
        integer :: cell_min_height(node%nreg)
        
        !*****
    
        !*****
        real(8) orrig_vol

        integer :: bottom = 0
        integer :: top = 0
        integer :: nextbottom = 0

		integer :: max_cl_height = 10

        integer :: error_int
        integer :: k, i, j
        
        integer :: nc
        integer :: start_height
        integer flag

		! this array will be used temperarily while enlarging the "collcell" array
        type(collcell_type), pointer, dimension(:) :: cccopy 
        !*****

        rmin = node%rmin
        rmax = node%rmax

        drfine = node%drfine
        dzfine = node%drfine
        nreg = node%nreg 

        do i = 1,nreg
            smallest_cc(i) = node%zregions(i)%smallest_cc
        end do


        do i = 1,nreg
            regzmin(i) = node%zregions(i)%zmin
            regzmax(i) = node%zregions(i)%zmax
        end do


        ! I'm leaving it to the user to make sure that the "smallest_cc" can be
        ! divided evenly into "reg_width".  If not, then the last
        ! cell's width will be extra thick, and that's the user's fault.

        ! this is the same for all regions
        reg_height = rmax - rmin


        ! working out widths in Z of the collision cells, for each zregion:

		! 'reg_width', 'cell_width', and 'num_wide' are all integers.


        do i = 1,nreg
			! "reg_width" is the number of dzfine units wide the region is.
            reg_width(i) = regzmax(i) - regzmin(i)

            ! the width, in dzfine units, of each collisoin cell in this zregion.
            cell_width(i) = int( smallest_cc(i) / dzfine )
			cell_width_real(i) = smallest_cc(i) / dzfine

        end do

        ! working out the minimum height in R (in drfine units) of the collision cells, for each zregion:
        ! maybe this won't ever be used.
        do i = 1,nreg
            cell_min_height(i) = even_round( smallest_cc(i) / drfine )
        end do

        ! initialize the "spec_r_marker" array
		do i = 1,nreg
        	if( node%zregions(i)%num_special_r > 0 ) then
            	allocate( node%zregions(i)%spec_r_marker(node%zregions(i)%num_special_r) )
				do j = 1,node%zregions(i)%num_special_r
            		node%zregions(i)%spec_r_marker(j) = 0
				end do
        	end if
		end do

        ! ALLOCATING THE COLLCELLS ARRAY
		node%ccarraysize = 1000
        allocate( node%collcells(node%ccarraysize), STAT=error_int )
        if( error_int .ne. 0 ) then
            write(*,*) "ERROR ALLOCATING ARRAY!  STOPPING PROGRAM"
            stop 928
        end if




		!*******************************
		!*******************************
		!*******************************



        ! coll cell counter
        nc = 0

        do k = 1,nreg

			! the target "volume"-like value which is aimed for in the creation of all collision cells.
            orrig_vol = dble(node%glob_rmax * drfine)**2 - (dble(node%glob_rmax * drfine) - smallest_cc(k))**2


			leftedge = regzmin(k) 
			rightedge = regzmin(k) + cell_width(k)
			width_remnant = cell_width_real(k) - dble(cell_width(k))

			if( rightedge > regzmax(k) .or. abs(dble(regzmax(k) - rightedge)) < cell_width_real(k) ) then
				rightedge = regzmax(k)
			end if

            do while( rightedge <= regzmax(k) .and. leftedge < regzmax(k) )

                bottom = rmax
                top = rmax 

        		! re-initialize the "spec_r_marker" array
				! we must do this every time we start again at the top of the column
        		if( node%zregions(k)%num_special_r > 0 ) then
					do j = 1,node%zregions(k)%num_special_r
            			node%zregions(k)%spec_r_marker(j) = 0
					end do
        		end if

                do while( bottom > rmin )

                    top = bottom
                    if( dble(top*drfine)**2-orrig_vol >= 0.0 ) then

                        bottom = even_round( sqrt(dble(top*drfine)**2-orrig_vol) / drfine )

                        if( (top - bottom) > max_cl_height ) then
                            bottom = top - max_cl_height
                        end if


                        if( bottom < rmin ) then
                            ! this happens at just the bottom of the processor
                            bottom = rmin
                        end if

                    else
                        ! this happens at the bottom of the entire simulation region
						if( top - rmin > max_cl_height ) then
							bottom = top - max_cl_height
						else
                        	bottom = rmin
						end if
                    end if



                    ! checking the next point to be calculated, and evening out the last couple of 
                    ! cells so we don't get a big and then a really really tiny cell.
                    if(  dble(bottom*drfine)**2-orrig_vol >= 0.0 ) then

                        nextbottom = even_round( sqrt(dble(bottom*drfine)**2-orrig_vol) / drfine )

                        if( (bottom - nextbottom) > max_cl_height ) then
                            nextbottom = bottom - max_cl_height
                        end if


                        if( nextbottom < rmin ) then
                            ! this happens at just the bottom of the processor
                            nextbottom = rmin
                        end if

                    else
                        ! this happens at the bottom of the entire simulation
						if( bottom - rmin > max_cl_height ) then
							nextbottom = bottom - max_cl_height
						else
                        	nextbottom = rmin
						end if
                    end if



                    ! adjust if the "nextbottom" is rmin, and it will end up making a tiny last cell
                    ! this is overridden if there's a "special_r" which readjusts the "bottom" value, below
                    if( nextbottom == rmin ) then
                        if( bottom > rmin      .and.   abs(bottom - rmin) < abs(top - bottom) &
														                    .and.   abs(top - rmin) > 14 ) then

                            bottom = abs(top - (abs(top - rmin) / 2))

                        else if( bottom > rmin .and.   abs(bottom - rmin) < abs(top - bottom) &
															                .and.   abs(top - nextbottom) <= 14 ) then

                            bottom = rmin

                        else
                            ! don't do anything; no problem.  nextbottom will be rmin is fine.
                        end if
                    end if



!					! preventing an extremely tiny ending
!					if( nextbottom == rmin ) then
!						bottom = rmin
!					end if




                    ! make sure there isn't a nearby "special_r" to adjust to.
                    if( node%zregions(k)%num_special_r > 0 ) then
                        do j = 1,node%zregions(k)%num_special_r

                            if( node%zregions(k)%special_r(j) == bottom ) then
                                node%zregions(k)%spec_r_marker(j) = 1 ! it's marked as having been lined up with.


                            else if(       node%zregions(k)%special_r(j) < top &
                                     .and. node%zregions(k)%special_r(j) > bottom &
                                     .and.  node%zregions(k)%spec_r_marker(j) == 0 ) then
                                ! a "special_r" is between top and bottom
                                bottom = node%zregions(k)%special_r(j)
                                node%zregions(k)%spec_r_marker(j) = 1 ! mark it

                            else if(       node%zregions(k)%special_r(j) < bottom &
                                     .and. node%zregions(k)%special_r(j) > nextbottom &
                                     .and. node%zregions(k)%spec_r_marker(j) == 0 ) then

                                ! a special_r is between bottom and nextbottom
                                ! we only shift "bottom" if it's closer to "bottom" than "nextbottom"
                                ! if it's closer to "nextbottom", then we'll wait until next time to shift (so we don't shift to make a tiny cell)
                                if( abs(bottom - node%zregions(k)%special_r(j)) <= abs(nextbottom - node%zregions(k)%special_r(j)) ) then
                                    bottom = node%zregions(k)%special_r(j)
                                    node%zregions(k)%spec_r_marker(j) = 1
                                end if

                            end if
                        end do
                    end if

				    !now we know the "top" and "bottom" values for this cell for cirtain
                    nc = nc + 1

                    if( nc > node%ccarraysize ) then
                        ! the array is full.  We must enlarge it.
					    allocate(cccopy(node%ccarraysize))
                        do j = 1,node%ccarraysize
                            cccopy(j) = node%collcells(j)
                        end do
                        deallocate(node%collcells)
                        allocate(node%collcells(node%ccarraysize*2))
                        do j = 1,node%ccarraysize
                            node%collcells(j) = cccopy(j)
                        end do
                        deallocate(cccopy)
                        node%ccarraysize = node%ccarraysize * 2
                    end if

                    node%collcells(nc)%rmin = bottom
                    node%collcells(nc)%rmax = top
                    node%collcells(nc)%zmin = leftedge
                    node%collcells(nc)%zmax = rightedge
                    ! this is where the extra wideness comes in if the smallest_cc is done wrong
                    if( (regzmax(k) - node%collcells(nc)%zmax) < cell_width(k) ) then
                        node%collcells(nc)%zmax = regzmax(k)
                    end if 

					! "psudoindex" is the line number you will find this cell in the "collcells.dat" file
					node%collcells(nc)%psudoindex = cellcounter_prev_sum + nc

					! load the volume, and "volcorrected" will be specified to "true" or "false", for whether it was corrected because
					! of metal interstection or not.

                    flag=0 ! collision cell volume correction
                    node%collcells(nc)%volume = corrected_cell_volume( dble( bottom ), &
                           dble( top  ), &
                           dble( node%collcells(nc)%zmin ), &
                           dble( node%collcells(nc)%zmax ), &
			   drfine, dzfine, node%polygons, node%npoly, node%collcells(nc)%volcorrected, node%collcells(nc)%psudoindex, flag)

					! for now, load a "0" for the mask
                    node%collcells(nc)%mask = 0



 
                end do

				! figureing out the cell width of the next iteration around

				leftedge = rightedge
				rightedge = rightedge + cell_width(k)
				width_remnant = width_remnant + (cell_width_real(k) - dble(cell_width(k)))
				if( int(width_remnant) > 0 ) then
					rightedge = rightedge + int(width_remnant)
					width_remnant = width_remnant - int(width_remnant)
				end if

				if( rightedge > regzmax(k) .or. abs(dble(regzmax(k) - rightedge)) < cell_width_real(k) ) then
					rightedge = regzmax(k)
				end if

				! just for a check
				if( leftedge == rightedge .and. rightedge < regzmax(k) ) then
					write(*,*) "ERROR! Somehow the next collcell's width will be zero"
					write(*,*) leftedge, "leftedge"
					write(*,*) rightedge, "rightedge"
					stop 1039
				end if


            end do

        end do


        node%nc = nc

        write(*,*) "Processor's number of collision cells:", node%nc


    end subroutine make_coll_cells












    !----even_round
    integer function even_round( num )
        implicit none
        real(8) num
        
        if( num - int(num) >= 0.5 ) then
            even_round = int(num) + 1
        else
            even_round = int(num)
        end if
 
    end function even_round










end module geombuilder_mod





