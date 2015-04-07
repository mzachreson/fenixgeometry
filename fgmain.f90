! FENIX GEMOETRY 10 -- to go with FENIX 10
! July, 2009, by Steven Schmidt and Ross Spencer
! Taken from previous versions of Fenix Geometry


! Total File list:

! fgtypes.f90
! fgdata.f90
! fgvolumecalc.f90
! fgbuilder.f90
! fgwrite.f90
! fgmain.f90

program fenixgeometry

    use types_mod
    use data_mod
    use write_mod
    use geombuilder_mod

    implicit none



    write(*,*) "Reading input_data.txt..."
    call read_input_data
    write(*,*) "Done reading input_data.txt"


    write(*,*) "Processing Geometry Data..."
    !call test_input
    call process_geom_data
    write(*,*) ""
    write(*,*) "Done processing geometry data"

    write(*,*) "Writing files..."
    call write_geom_files
    write(*,*) "Done writing files."


end program fenixgeometry



