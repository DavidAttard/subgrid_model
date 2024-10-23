program cond_mass_funct_david

  ! cond_mass_funct.f90 - serially reads a distributed checkpoint and
  ! halo catalogue and produces a conditional mass function (halo
  ! counts binned in 10 mass bins) and collapsed fraction vs. overdensity 
  ! in equal-volume cubical regions of a given size.
  !
  ! Compile with:
  ! ifort -fpp -DBINARY -mcmodel=medium -shared-intel cond_mass_funct_david.f90 -o cond_mass_funct_david -O3
  !
  ! by Ilian Iliev : Oct., 2009 

  use precision  
  implicit none
  

  !list of redshift to generate slices for 

  character(len=*),parameter :: checkpoints='./../checkpoints_halos' 

  !cosmology
  real*8,parameter :: omegabh2=0.0219961,omega0=0.318,lambda0=0.682,h=0.678&
       & ,an=0.9611,tcmb=2.726,s8=0.833 

  ! frequently used parameters are found in this header file:
  include '/research/prace5/shared/iti20/IDM/cubep3m_200213_8_4096_100Mpc_ext2_pencil/workdir/parameters'

  real*8 :: mass0, m_grid, m_min, m_max, m_box
  real*8 :: rho_crit_0, rho_crit_z, rho_crit_0_MpcM_sun, rho_bar, rho_bar_cgs

  real*8,parameter :: M_sun = 1.989d33, cm_to_Mpc=3.086d24

  ! Halo mass binning

  real(kind=4), parameter :: lM_min = 9  ! (log10) minimum halo
  ! mass for binning
  real(kind=4), parameter :: lM_max = 9.4 ! (log10) maximum halo
  ! mass for binning. Just a single bin for all halos above M_max
  real(kind=4), parameter :: n_bin = 2     ! total number of mass
  ! bins; 0.5= mass bin size (in log10)
  real(kind=4) :: mms(n_bin+1)  

  ! size of (cubical) regions to calculate the binned mass function
  ! and density in. Units are internal code units (fine grid cells).

  integer(kind=4), parameter :: n_region = 256  ! number of coarse
!  integer(kind=4), parameter :: n_region = 256  ! number of coarse
  ! regions per dimension. This value corresponds to 1 RT cell of 
  ! 114/h Mpc box at 256^3 cells and box=114/h N-body
  ! box to be analyzed (last value is supplied by 'parameters')


  real(kind=8) :: coarse_halos(n_region,n_region,n_region,n_bin+1,0:1) 
  ! this array will store the coarsened halo field
  ! last index is 0 for n_halos and 1 for m_halos 
  ! Note that this array will become huge if
  ! n_region is a large number!!!
  real(kind=8) :: coll_frac(n_region,n_region,n_region)  ! this array
  ! will store the local collapsed fraction

  !Halo data arrays                                                  
  !      
  real(4), dimension(3) :: halo_pos ! halo position (cells)
  real(4), dimension(3) :: x_mean ! centre of mass position 
  real(4), dimension(3) :: v_mean ! velocity
  real(4), dimension(3) :: l      ! angular momentum
  real(4) :: v_disp      ! velocity dispersion
  real(4) :: radius_calc ! radius
  real(4) :: halo_mass,halo_mass_uncorrected ! mass calculated on the
  ! grid (in grid masses)
  real(4) :: imass       ! mass calculated from the particles (in
  ! grid masses)
  integer(kind=4) :: dim 


  real(kind=4), parameter :: ncc = nc/nodes_dim 
  !linear size of 1 sub-volume, in cells 

  integer,dimension(3) :: p

  integer(kind=4) :: counter = 0 

  integer(kind=4), parameter :: nn = nodes_dim**3 !total number of
  ! sub-volumes

  ! the rest are all internal variables

  character(len=512) :: ifile1, ifile2, ofile1, ofile2, ofile3, ofile4
  character(len=7) :: z_s
  character(len=4) :: rank_s
  character(len=4) :: size_s
  integer(4), parameter :: max_input=100

  real(4) :: a,t,tau,dt_f_acc,dt_c_acc,mass_p,dt_pp_acc

  integer(4) :: np_local,np_current, nts, sim_checkpoint,&
       & sim_projection, sim_halofind

  integer(4), parameter :: max_np = density_buffer * ( ((nf_tile-2 &
       &*nf_buf)*tiles_node_dim/2)**3 + (8*nf_buf**3 + 6*nf_buf &
       &*(((nf_tile-2*nf_buf)*tiles_node_dim)**2) + 12*(nf_buf**2) &
       &*((nf_tile-2*nf_buf)*tiles_node_dim))/8.0 )
  
  
  integer, parameter :: nhv=23  !! number of halo quantities in each (NOTE there are actually 22 but in the text file each entry takes up
  !! two lines with the 2nd line being indented. This indentation is marked as a 0 when reading and thus to read all of the entries you need to
  !! read one more character)
  !! catalog entry
  integer*4 :: n_halo
  real(4), dimension(nhv) :: halo_input_buffer

  integer(kind=4) :: rank,i,j,k,cp,nploc,fstat,num_checkpoints, ii
  real(kind=4) :: z
  real(kind=4), dimension(max_input) :: z_checkpoint
  integer :: node_coords(nn,3)           !! coordinates of each node
  !! in global cube

  ! density in file is in real(kind=4), read in via this array
  real :: ndens_real(n_region,n_region,n_region)
  real(kind=dp) :: coarse_dens(n_region,n_region,n_region)

  character(len=180) :: dens_file
  real(kind=dp) :: avg_dens

  rho_crit_0=1.88d-29*h**2!g/cm^3
  rho_crit_0_MpcM_sun=rho_crit_0*cm_to_Mpc**3/M_sun
  rho_bar = rho_crit_0_MpcM_sun*omega0
  rho_bar_cgs = rho_crit_0*omega0
  print*,'rho_bar(Mpc^3/Msun)=',rho_bar
  print*,'rho_bar_cgs=',rho_bar_cgs
  !
  M_box = rho_bar*(box/h)**3
  M_grid = M_box/real(nc)**3

  !set up mass bins

  m_min = 10**lM_min

  do i=1,n_bin+1
     mms(i)=m_min*10.**((lM_max-lM_min)*(real(i-1))/n_bin)!mass bins to sort halos
     ! into
  end do

  !! Calculate node_coords
  do rank=0,nn-1
     do k=1,nodes_dim
        do j=1,nodes_dim
           do i=1,nodes_dim
              if (rank == (i-1)+(j-1)*nodes_dim+(k-1)*nodes_dim**2)&
                   & node_coords(rank,:)=(/(i-1),(j-1),(k-1)/)
           enddo
        enddo
     enddo
  enddo
  print*,'Max node_coords',node_coords(nn-1,1:3)
  print*,'nc ',nc
  print*,'check ncc',ncc

  !! Read in list of checkpoints to create the data for

  open(11,file=checkpoints,status='old',iostat=fstat)
  if (fstat /= 0) write(*,*) 'error opening checkpoint list file:' &
       &,checkpoints
  do num_checkpoints=1,max_input
     read(unit=11,err=51,end=41,fmt='(f20.10)')&
          & z_checkpoint(num_checkpoints)
  enddo
41 num_checkpoints=num_checkpoints-1
51 close(11)
  write(*,*) 'checkpoints to slice:'
  do i=1,num_checkpoints
     write(*,'(f5.1)') z_checkpoint(i)
  enddo


  write(size_s,'(i4)') n_region
  size_s=adjustl(size_s)  !this will label the coarse grid size in
  ! filenames

  !Loop over checkpoints

  do cp=1,num_checkpoints
     write(z_s,'(f7.3)') z_checkpoint(cp)
     z_s=adjustl(z_s)

     ofile1='./../binning/'//z_s(1:len_trim(z_s)) &
          &//'dens'//size_s(1:len_trim(size_s))//'.bin'
     ofile2='./../binning/'//z_s(1:len_trim(z_s)) &
          &//'halo_num'//size_s(1:len_trim(size_s))//'.bin'
     ofile3='./../binning/'//z_s(1:len_trim(z_s)) &
          &//'halo_mass'//size_s(1:len_trim(size_s))//'.bin'
     ofile4='./../binning/'//z_s(1:len_trim(z_s)) &
          &//'coll_frac'//size_s(1:len_trim(size_s))//'.bin'


     !! open output_slice

     open(unit=12,file=ofile1,form='unformatted')
     open(unit=13,file=ofile2)
     open(unit=16,file=ofile4,form='unformatted')
     open(unit=17,file=ofile3)

     coarse_halos=0.0d0
     coll_frac=0.0d0

     !! loop over each rank 
      print*,nn
      
     do rank=0,0
        write(rank_s,'(i4)') rank
        rank_s=adjustl(rank_s)
        !important to use the ntot_all (halos not exluded) file and not the n_all (halos_excluded) file 
        ifile1='/research/prace5/shared/iti20/IDM/sph_smooth_cubep3m_200213_8_4096_100Mpc_ext2_pencil/global/so/nc256/'//z_s(1:len_trim(z_s))//'ntot_all.dat'
        ifile2='/research/prace5/shared/iti20/IDM/cubep3m_200213_8_4096_100Mpc_ext2_pencil/results/full_halo_catalogues/p/scratch/chpo22/hpo222/cubep3m_200213_8_4096_100Mpc_ext2_pencil/results/'//z_s(1:len_trim(z_s))//'halo.dat'

       if (fstat /= 0) then
           write(*,*) 'error opening catalog'
           write(*,*) 'rank',rank,'file:',ifile2
           stop !call mpi_abort(mpi_comm_world,ierr,ierr)
        endif


        open(unit=15,file=ifile2, form='formatted')!full halo catalogues (ascii)

         ! Open density file: note that it is in `unformatted' form
         open(unit=14,file=ifile1,form="binary", access="sequential",status="old")  

         read(14) dim,dim,dim 
         ! Read in data and store it in ndens
         read(14) ndens_real
         coarse_dens(:,:,:)=ndens_real(:,:,:)
               
        !we are done with the particles for this checkpoint
        close(14)

        print*,'Total mass:' ,sum(coarse_dens)
        print*,'maximum density: ',maxval(coarse_dens)
        print*,minval(coarse_dens)
        !now do the halos

        !read halo file header

        !read(15) n_halo !distributed halo catalogues only
        !print*, n_halo

         !print*,'check halo data'
        do 
           read(15,end=112,err=113, fmt='(22f20.10)') halo_input_buffer
                !halo_pos(:), x_mean(:), v_mean(:), l(:), v_disp,&
                !radius_calc, halo_mass, imass, halo_mass_uncorrected
                !print*, halo_input_buffer
                

                halo_pos = halo_input_buffer(1:3)
                halo_mass = halo_input_buffer(4) !the defintiions are found in halofind_particles.f90 line 305
                
                
                      !print*,halo_pos(:), x_mean(:), v_mean(:), l(:), v_disp,&
                       !    radius_calc, halo_mass, imass, halo_mass_uncorrected
           
                     
           !           halo_input_buffer(1:3)=modulo(halo_input_buffer(1:3),ncc)

           !put halo centers in local coords
           !halo_pos(1:3)=modulo(halo_pos(1:3),ncc)

           !now set them back to global coordinates, but now ensuring
           !they are the same as the ones for the particles 

           !halo_pos(1:3) = halo_pos(1:3)+node_coords(rank,1:3)*ncc

           !position of halo center in the coarse grid (we ignore
           !here that parts of the halo that might be in another cuboid, 
           !i.e. all halo mass is assigned to the coarse cell containing 
           !the halo centre) 

           p(1:3)=modulo(int(halo_pos(1:3)*real(n_region)/real(nc)),n_region)+1
           
           mass0=halo_mass*M_grid           

           if((mass0<=10**lM_min).or.(mass0>10**lM_max))then
              !print*,'Halo mass not within bounds'
           else
              if(mass0>=mms(n_bin+1))then !halo is large, above 10^12
                 ! M_solar
                 coarse_halos(p(1),p(2),p(3),n_bin+1,0) &
                      &=coarse_halos(p(1),p(2),p(3),n_bin+1,0)+1
                 coarse_halos(p(1),p(2),p(3),n_bin+1,1) &
                      &=coarse_halos(p(1),p(2),p(3),n_bin+1,1)+mass0
              else !find which mass bin the halo falls in
                 do ii=2,n_bin+1 
                    if(mass0>mms(ii-1).and.mass0<=mms(ii)) then
                       coarse_halos(p(1),p(2),p(3),ii-1,0) &
                            &=coarse_halos(p(1),p(2),p(3),ii-1,0)+1
                       coarse_halos(p(1),p(2),p(3),ii-1,1) &
                            &=coarse_halos(p(1),p(2),p(3),ii-1,1)+mass0
                       !print*,'check 2 ',mass0,mms(ii),mms(ii+1),ii
                    end if
                  end do
            end if
            coll_frac(p(1),p(2),p(3))=coll_frac(p(1),p(2),p(3)) + mass0
           end if
        enddo
113        print*,' Halo catalogue ',z_s(1:len_trim(z_s)),rank,' gives a reading error.'     

112     close(15)

     end do
     
     !fix units - go to CGS

     coarse_dens=coarse_dens*M_grid
     do i=1,n_region
        do j=1,n_region
           do k=1,n_region
              if(coarse_dens(i,j,k)>0)then
                 coll_frac(i,j,k)=coll_frac(i,j,k)/coarse_dens(i,j,k)
              else
                 coll_frac(i,j,k)=0.0
              end if
           end do
        end do
     end do

     !Density diagnostics

     print*,'Total mass:' ,sum(coarse_dens)
     print*,'Total mass, theory: ',M_box
     print*,'minimum density: ',minval(coarse_dens)
     print*,'maximum density: ',maxval(coarse_dens)
     print*,'average density: ',sum(coarse_dens)/n_region**3

     print*,'Total coll. fraction:' ,sum(coll_frac*coarse_dens)/sum(coarse_dens),z_checkpoint(cp)
     print*,'minimum coll. fraction: ',minval(coll_frac)
     print*,'maximum coll. fraction: ',maxval(coll_frac)
     print*,'total number of halos: ',sum(coarse_halos(:,:,:,:,0))

     write(12) n_region,n_region,n_region
     write(12) (((real(coarse_dens(i,j,k),4),i=1,n_region),j=1,n_region),k=1,n_region)

     do i=1,n_region
        do j=1,n_region
           do k=1,n_region
              write(13,*) real(coarse_dens(i,j,k),4),(real(coarse_halos(i,j,k,ii,0),4),ii=1,n_bin+1)
              write(17,*) real(coarse_dens(i,j,k),4),(real(coarse_halos(i,j,k,ii,1),4),ii=1,n_bin+1),real(coll_frac(i,j,k),4)
           end do
        end do
     end do

     write(16) n_region,n_region,n_region
     write(16) (((real(coll_frac(i,j,k),4),i=1,n_region),j=1,n_region),k=1,n_region)

     write(*,*) 'done z=',z_checkpoint(cp)

  enddo
  
end program cond_mass_funct_david
