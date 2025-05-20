program cond_mass_funct_1024_256_mpi 

  ! cond_mass_funct.f90 - serially reads a distributed checkpoint and
  ! halo catalogue and produces a conditional mass function (halo
  ! counts binned in 10 mass bins) and collapsed fraction vs. overdensity 
  ! in equal-volume cubical regions of a given size.
  !
  ! Compile with:
  ! ifx -fpp -DBINARY -mcmodel=medium -shared-intel cond_mass_funct_1024_256_mpi.f90 -o cond_mass_funct_1024_256_mpi -O3
  !

  use precision  
  implicit none

  character(len=*),parameter :: checkpoints='./../checkpoints_halos'

  ! cosmology and parameters
  real*8,parameter :: omegabh2=0.0219961,omega0=0.318,lambda0=0.682,h=0.678&
       & ,an=0.9611,tcmb=2.726,s8=0.833 
  include './../parameters'

  real*8 :: mass0, m_grid, m_min, m_max, m_box
  real*8 :: rho_crit_0, rho_crit_z, rho_crit_0_MpcM_sun, rho_bar, rho_bar_cgs

  real*8,parameter :: M_sun = 1.989d33, cm_to_Mpc=3.086d24
  character(len=10) :: cg_s
  real(kind=4), parameter :: lM_min = 8, lM_max = 9, n_bin = 2
  real(kind=4) :: mms(n_bin+1)  

  integer(kind=4), parameter :: n_region = 1024, ng = 1024
  integer, parameter :: coarsened_grid = 256, block = ng / coarsened_grid
  integer :: i, j, k, ii, jj, kk, dx, dy, dz, ci, cj, ck, shift_id, fh
  real(kind=dp) :: coarse_dens_coarse(n_region,n_region,n_region), coll_frac_coarse(n_region,n_region,n_region)
  character(len=50) :: fname
  real(kind=8) :: coll_frac(n_region,n_region,n_region)
  real(4), dimension(3) :: halo_pos
  real(4) :: halo_mass, halo_mass_uncorrected, imass
  integer(kind=4) :: dim
  real(kind=4), parameter :: ncc = nc/nodes_dim 
  integer,dimension(3) :: p
  integer(kind=4), parameter :: nn = nodes_dim**3
  character(len=512) :: ifile1, ifile2, ofile1, ofile2, ofile3, ofile4
  character(len=7) :: z_s
  character(len=4) :: rank_s, size_s
  integer(4), parameter :: max_input=100
  real(4) :: a,t,tau,dt_f_acc,dt_c_acc,mass_p,dt_pp_acc
  integer(4) :: np_local,np_current, nts, sim_checkpoint, sim_projection, sim_halofind
  integer(4), parameter :: max_np = density_buffer * ( ((nf_tile-2 *nf_buf)*tiles_node_dim/2)**3 + (8*nf_buf**3 + 6*nf_buf *( (nf_tile-2*nf_buf)*tiles_node_dim)**2 + 12*(nf_buf**2) *(nf_tile-2*nf_buf)*tiles_node_dim)/8.0 )
  integer, parameter :: nhv=23
  integer*4 :: n_halo
  real(4), dimension(nhv) :: halo_input_buffer
  integer(kind=4) :: rank,cp,nploc,fstat,num_checkpoints
  real(kind=4) :: z
  real(kind=4), dimension(max_input) :: z_checkpoint
  integer :: node_coords(nn,3)
  real :: ndens_real(n_region,n_region,n_region)
  real(kind=dp) :: coarse_dens(n_region,n_region,n_region)
  character(len=180) :: dens_file
  real(kind=dp) :: avg_dens
  integer :: cp_index

character(len=32) :: arg
integer :: ios

call get_command_argument(1, arg)
read(arg, *, iostat=ios) cp_index
cp_index = cp_index
if (ios /= 0) then
    print *, 'Error: invalid or missing checkpoint index argument.'
    stop
endif

  !--- COSMOLOGY CALCULATIONS ---
  rho_crit_0=1.88d-29*h**2
  rho_crit_0_MpcM_sun=rho_crit_0*cm_to_Mpc**3/M_sun
  rho_bar = rho_crit_0_MpcM_sun*omega0
  rho_bar_cgs = rho_crit_0*omega0
  print*,'rho_bar(Mpc^3/Msun)=',rho_bar
  print*,'rho_bar_cgs=',rho_bar_cgs
  M_box = rho_bar*(box/h)**3
  M_grid = M_box/real(nc)**3

  m_min = 10**lM_min
  do i=1,n_bin+1
     mms(i)=m_min*10.**((lM_max-lM_min)*(real(i-1))/n_bin)
  end do

  ! Calculate node_coords
  do rank=0,nn-1
     do k=1,nodes_dim
        do j=1,nodes_dim
           do i=1,nodes_dim
              if (rank == (i-1)+(j-1)*nodes_dim+(k-1)*nodes_dim**2) &
                 node_coords(rank,:)=(/(i-1),(j-1),(k-1)/)
           enddo
        enddo
     enddo
  enddo

  ! Read list of redshift checkpoints
  open(11,file=checkpoints,status='old',iostat=fstat)
  if (fstat /= 0) stop 'error opening checkpoint list file'
  do num_checkpoints=1,max_input
     read(unit=11,err=51,end=41,fmt='(f20.10)') z_checkpoint(num_checkpoints)
  enddo
41 num_checkpoints=num_checkpoints-1
51 close(11)

  if (cp_index < 1 .or. cp_index > num_checkpoints) then
     print *, 'Invalid checkpoint index. Exiting.'
     stop
  end if

  ! Proceed with selected checkpoint
  cp = cp_index
  write(*,*) 'Processing z=',z_checkpoint(cp)
  write(z_s,'(f7.3)') z_checkpoint(cp)
  z_s=adjustl(z_s)
  write(size_s,'(i4)') n_region
  write(cg_s, '(I0)') coarsened_grid
  size_s=adjustl(size_s)


     write(*,*) 'doing z=',z_checkpoint(cp)
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


     coll_frac=0.0d0

     !! loop over each rank 
      print*,nn
      
     do rank=0,0
        write(rank_s,'(i4)') rank
        rank_s=adjustl(rank_s)
        !important to use the ntot_all (halos not exluded) file and not the n_all (halos_excluded) file 
        ifile1='/p/largedata2/hpo22/subgrid_modeling/sph_smooth_cubep3m_200213_8_4096_100Mpc_ext2_pencil/global/so/nc1024/'//z_s(1:len_trim(z_s))//'ntot_all.dat'
        ifile2 = '/p/largedata2/hpo22/subgrid_modeling/cubep3m_200213_8_4096_100Mpc_ext2_pencil/results/' & 
       // 'full_halo_catalogues/p/scratch/chpo22/hpo222/cubep3m_200213_8_4096_100Mpc_ext2_pencil/results/' & 
       // z_s(1:len_trim(z_s)) // 'halo.dat'
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
            coll_frac(p(1),p(2),p(3))=coll_frac(p(1),p(2),p(3)) + mass0
           end if
        enddo
113        print*,' Halo catalogue ',z_s(1:len_trim(z_s)),rank,' gives a reading error.'     

112     close(15)

     end do
     
     !fix units - go to CGS

     coarse_dens=coarse_dens*M_grid

     print*, 'First coarse dens value: ', coarse_dens(1,1,1)
     print*, 'First coll frac value: ', coll_frac(1,1,1)
      
      shift_id = 0
      do dz = 0, block-1
         do dy = 0, block-1
            do dx = 0, block-1
            shift_id = shift_id + 1
            write(fname, '(A,I2.2,A,A,A)') './../binning/coarsened_shift_', shift_id, '_', trim(cg_s)//'_z'//trim(z_s), '.dat'

            ! Initialize coarse arrays
            coarse_dens_coarse = 0.0
            coll_frac_coarse = 0.0

            ! Coarsen: loop over coarse grid cells
            do ci = 1, coarsened_grid
               do cj = 1, coarsened_grid
                  do ck = 1, coarsened_grid
                  do ii = 0, block-1
                     do jj = 0, block-1
                        do kk = 0, block-1
                        i = modulo((ci-1)*block + ii + dx, ng) + 1
                        j = modulo((cj-1)*block + jj + dy, ng) + 1
                        k = modulo((ck-1)*block + kk + dz, ng) + 1

                        coarse_dens_coarse(ci,cj,ck) = coarse_dens_coarse(ci,cj,ck) + coarse_dens(i,j,k)
                        coll_frac_coarse(ci,cj,ck) = coll_frac_coarse(ci,cj,ck) + coll_frac(i,j,k)
                        end do
                     end do
                  end do

                  end do
               end do
            end do

            print *, ' Saving the coarsened grid: coarse_dens_coarse(1,1,1)', coarse_dens_coarse(1,1,1)
            print *, ' Saving the coarsened grid: coll_frac_coarse(1,1,1)', coll_frac_coarse(1,1,1)
            ! Write to file
            open(unit=fh, file=fname, form='formatted')
            do ci = 1, coarsened_grid
               do cj = 1, coarsened_grid
                  do ck = 1, coarsened_grid
                  write(fh, *)  real(coarse_dens_coarse(ci,cj,ck)), real(coll_frac_coarse(ci,cj,ck))
                  end do
               end do
            end do
            close(fh)
            print *, 'Wrote: ', trim(fname)

            end do
         end do
      end do

     write(*,*) 'done z=',z_checkpoint(cp)

end program cond_mass_funct_1024_256_mpi
