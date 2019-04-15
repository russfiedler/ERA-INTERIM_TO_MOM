program make_era_interim
!
! Usage::
! make_era_interim -v var -x var -i file_in -o file_out -a [yYnN] -f [yYnN]
!
! -v var      : name of variable to write 
! -x var      : Variable to extract from original
! -i file_in  : input file
! -o file_out : output file
! -a y/n      : is input variable accumulated or not
! -f y/n      : do we want to mask out land values and fill with nearby ocean values.
! 
use iso_fortran_env
use netcdf
implicit none
integer(int16),allocatable,dimension(:,:,:) :: var_in
real(real32),allocatable,dimension(:,:,:)   :: var_out
real(real64),allocatable,dimension(:)   :: time
integer(int32),allocatable,dimension(:)   :: timein
integer(int32),allocatable,dimension(:,:)   :: mask
real(real64) :: scale_factor,add_offset
integer(int32) :: numargs
character*132 :: file_in,file_out
character*1   :: fill, accum
character*2   :: flag
character*32 cvar,cvarin,timename
integer :: i,j,k
integer :: ncidin,tidin,vidin,did,dimlen,vidout,ncidout,tidout

numargs = command_argument_count()
if ( numargs /= 12 ) then
   write(*,*) 'Usage: make_era_interim -v var -x varin -i file_in -o file_out -f [y/n] -a accum [y/n]'
   stop 1
endif

do i = 1,11,2
   call get_command_argument(i,flag)
   select case(trim(flag))
   case('-v')
      call get_command_argument(i+1,cvar)
   case('-x')
      call get_command_argument(i+1,cvarin)
   case('-i')
      call get_command_argument(i+1,file_in)
   case('-o')
      call get_command_argument(i+1,file_out)
   case('-f')
      call get_command_argument(i+1,fill)
   case('-a')
      call get_command_argument(i+1,accum)
   case DEFAULT
      write(*,*) 'Usage: make_era_interim -v var -x varin -y year -i file_in -o file_out -f [y/n] -a accumulated [y/n]'
      stop 1
   end select
enddo

call handle_netcdf(nf90_open(trim(file_in),nf90_nowrite,ncidin))
call handle_netcdf(nf90_inq_dimid(ncidin,'time',did))
call handle_netcdf(nf90_inquire_dimension(ncidin,did,len=dimlen))
call handle_netcdf(nf90_inquire_dimension(ncidin,did,name=timename))
call handle_netcdf(nf90_inq_varid(ncidin,trim(timename),tidin))
call handle_netcdf(nf90_inq_varid(ncidin,trim(cvarin),vidin))
write(*,*)  'Got input files'
call handle_netcdf(nf90_open(trim(file_out),nf90_write,ncidout))
call handle_netcdf(nf90_inq_varid(ncidout,'Time',tidout))
call handle_netcdf(nf90_inq_varid(ncidout,trim(cvar),vidout))
write(*,*)  'Got output files'


! Convert time to days
allocate(timein(dimlen),time(dimlen))
allocate(var_out(240,121,dimlen))
allocate(var_in(240,121,dimlen))
call handle_netcdf(nf90_get_var(ncidin,tidin,timein))
time = timein/24.0
call handle_netcdf(nf90_put_var(ncidout,tidout,time))
call handle_netcdf(nf90_get_var(ncidin,vidin,var_in))
call handle_netcdf(nf90_get_att(ncidin,vidin,'scale_factor',scale_factor))
call handle_netcdf(nf90_get_att(ncidin,vidin,'add_offset',add_offset))
call handle_netcdf(nf90_close(ncidin))

print *, 'Read data'

! Idiotic reversal of axes. Sigh.
do k = 1,dimlen
   do j=1,121
      var_out(:,j,k) = real(scale_factor,real32)*var_in(:,122-j,k)+real(add_offset,real32)
   enddo
enddo
deallocate(var_in)
print *, 'Scaled data'

! if accumulated we actually want v1],v[2]-v[1],v[3-v2],v[4]-v[3]
if (accum == 'y' .or. accum == 'Y' ) then
   do k = 1,dimlen,4
      var_out(:,:,k+3)=(var_out(:,:,k+3)-var_out(:,:,k+2))/(3600.0*3)
      var_out(:,:,k+2)=(var_out(:,:,k+2)-var_out(:,:,k+1))/(3600.0*3)
      var_out(:,:,k+1)=(var_out(:,:,k+1)-var_out(:,:,k))/(3600.0*3)
      var_out(:,:,k)=var_out(:,:,k)/(3600.0*3)
   enddo
endif

if ( fill == 'y' .or. fill == 'Y' ) then
   write(*,*) 'Filling land'
   call handle_netcdf(nf90_open('era_land_sea_extra.nc',nf90_nowrite,ncidin))
   call handle_netcdf(nf90_inq_varid(ncidin,'lsm',vidin))
   allocate(mask(240,121))
   call handle_netcdf(nf90_get_var(ncidin,vidin,mask))
   call handle_netcdf(nf90_close(ncidin))
   
   call fillerup(var_out,mask)
   deallocate(mask)
endif
      
   

print *, 'Writing'

call handle_netcdf(nf90_put_var(ncidout,vidout,var_out))

! End of interval
call handle_netcdf(nf90_close(ncidout))

contains

subroutine fillerup(var,mask)
implicit none
real,dimension(:,:,:),intent(inout) :: var
integer,dimension(:,:),intent(in) :: mask
real,allocatable,dimension(:,:) :: tmp,mask_new,mask_newer
integer :: i,j,k
integer :: i_in,j_in,its
integer :: j_lo,j_hi
real :: sm

i_in = size(mask,dim=1)
j_in = size(mask,dim=2)

allocate(tmp(i_in,j_in),mask_new(i_in,j_in),mask_newer(i_in,j_in))

do k=1,size(var,dim=3)
   print *,k,i_in,j_in,size(var,dim=3)
   mask_new = 1.-mask
   mask_newer=mask_new
   its = 0
   do 
      its = its+1
      tmp=var(:,:,k)*mask_newer
      mask_new = mask_newer
      do j = 1,j_in
         j_lo=max(j-1,1)
         j_hi=min(j+1,j_in)
         do i = 2,i_in-1
            if (mask_new(i,j) == 0.0 ) then
               sm = sum(mask_new(i-1:i+1,j_lo:j_hi))
               if(sm .ne. 0.0 ) then
                  var(i,j,k)= sum(tmp(i-1:i+1,j_lo:j_hi))/sm
                  mask_newer(i,j) = 1.0
               endif
            endif
         enddo
         i=1
         if (mask_new(i,j) == 0.0 ) then
            sm = sum(mask_new(i:i+1,j_lo:j_hi))+sum(mask_new(i_in,j_lo:j_hi))
            if(sm .ne. 0.0 ) then
               var(i,j,k)= (sum(tmp(i:i+1,j_lo:j_hi))+sum(tmp(i_in,j_lo:j_hi)))/sm
               mask_newer(i,j) = 1.0
            endif
         endif
         i=i_in
         if (mask_new(i,j) == 0.0 ) then
            sm = sum(mask_new(1,j_lo:j_hi))+sum(mask_new(i-1:i,j_lo:j_hi))
            if(sm .ne. 0.0 ) then
               var(i,j,k)= (sum(tmp(1,j_lo:j_hi))+sum(tmp(i-1:i,j_lo:j_hi)))/sm
               mask_newer(i,j) = 1.0
            endif
         endif
      enddo
      if(all(mask_newer == 1.0 ) .or. its .eq. 30 ) exit
   enddo
enddo
end subroutine fillerup

subroutine handle_netcdf(error_code,fatal)
integer, intent(in) :: error_code
logical, intent(in), optional :: fatal
! Placeholder
!integer, intent(in), dimension(:), optional :: nonfatal_code

logical :: isfatal

isfatal = .true.
if(present(fatal) ) then
  isfatal =fatal
endif
if(error_code /= nf90_noerr ) then
   write(*,*) 'Error: ',nf90_strerror(error_code)
   if(isfatal) stop
endif
end subroutine handle_netcdf


end program make_era_interim
         

         


