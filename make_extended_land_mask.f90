program make_extended_land_mask
! Make some sea point which are polluted by land values land
!
! Usage::
! make_extended_land_mask
! Use this on the standard ERA-INTERIM land mask file. I've renamed it.
! 
use iso_fortran_env
use netcdf
implicit none
integer(int32),allocatable,dimension(:,:)   :: mask
integer(int32),dimension(17) :: xind,yind
integer :: i
integer :: ncidin,vidin

! Points to remove
xind = (/ 137,111,112, 83, 82,134, 79, 81, 81,152,153,154,153,154, 96,120,108/)
yind = (/  74, 47, 46, 70, 70, 98, 55, 55, 54, 97, 97, 97, 96, 96, 92, 48, 55/)
call handle_netcdf(nf90_open('era_land_sea_extra.nc',nf90_write,ncidin))
call handle_netcdf(nf90_inq_varid(ncidin,'lsm',vidin))
write(*,*)  'Got input files'


allocate(mask(240,121))
call handle_netcdf(nf90_get_var(ncidin,vidin,mask))
do i=1,17
   mask(xind(i),yind(i))=1
enddo

print *, 'Writing'

call handle_netcdf(nf90_put_var(ncidin,vidin,mask))

! End of interval
call handle_netcdf(nf90_close(ncidin))

contains


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
!if(error_code /= nf90_noerror ) then
if(error_code /= nf90_noerr ) then
   write(*,*) 'Error: ',nf90_strerror(error_code)
   if(isfatal) stop
endif
end subroutine handle_netcdf


end program make_extended_land_mask
         

         


