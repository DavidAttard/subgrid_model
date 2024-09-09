program read_text_file
    implicit none

    integer :: ios
    real(8) :: halo_input_buffer(23)  ! Adjusted to hold 22 floating-point numbers
    character(len=512) :: ifile2
    real(8) :: max_value
    logical :: first_read

    ! Specify the file path
    ifile2 = '/research/prace5/shared/iti20/IDM/cubep3m_200213_8_4096_100Mpc_ext2_pencil/results/full_halo_catalogues/p/scratch/chpo22/hpo222/cubep3m_200213_8_4096_100Mpc_ext2_pencil/results/27.900halo.dat'
    ifile2='/research/prace5/shared/iti20/IDM/cubep3m_200213_8_4096_100Mpc_ext2_pencil/results/full_halo_catalogues/p/scratch/chpo22/hpo222/cubep3m_200213_8_4096_100Mpc_ext2_pencil/results/26.120halo.dat'
    ifile2='/research/prace5/shared/iti20/IDM/cubep3m_200213_8_4096_100Mpc_ext2_pencil/results/full_halo_catalogues/p/scratch/chpo22/hpo222/cubep3m_200213_8_4096_100Mpc_ext2_pencil/results/7.570halo.dat'
    ! Open the file for reading with unit number 15
    open(unit=15, file=ifile2, form='formatted', status='old')

    ! Begin the indefinite loop
    print*,modulo(int(8197*real(256)/real(8192)),256)+1
    print*,modulo(int(1*real(256)/real(8192)),256)+1
    print*,modulo(int(8000*real(256)/real(8192)),256)+1

    print*,modulo(8000,1024)
    print*,modulo(8000,1024)+7*1024

    do 
        ! Attempt to read a line of data into halo_input_buffer
        read(15,end=112,err=113,fmt='(22f20.10)') halo_input_buffer
        ! Update max_value if halo_input_buffer(1) is larger
        if (first_read .or. halo_input_buffer(1) > max_value) then
            max_value = halo_input_buffer(1)
            first_read = .false.
        end if
 
    enddo
    113        print*,' Halo catalogue gives a reading error.'     

    112     close(15)

    print*,'Max value: ',max_value
end program read_text_file


        !print*,'Max x_posn=', maxval(halo_input_buffer(1)) 
        
        !print*,'Min x_posn=', minval(halo_input_buffer(1))       
