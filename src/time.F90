!!****f* source/util/tools/fdate
!!
!! NAME
!!   ftime - formated date
!!
!! SYNOPSIS
!!   d = fdate()
!!   character (len=*) = fdate()
!!
!! DESCRIPTION
!!   function to return out the date, formated as.
!!		 "MM-DD-YY  HH:mm.ss".
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!***

module timers
implicit none



contains



function f_date_func ()

  character(len=24) :: f_date_func
  character(len=8)  :: date
  character(len=10) :: time
  character(len=5)  :: timezone
  integer           :: values(8)

  call date_and_time (date, time, timezone, values)

  write (f_date_func, &
       &   "(I2.2,'-',I2.2,'-',I4.4,' ',I2.2,':',I2.2,'.',I2.2,'     ')") &
       &   values(2), values(3), values(1), values(5), values(6), &
       &   values(7)

  return
end function f_date_func



subroutine time_message(unit,mesg)
  implicit none
  integer, intent(in) :: unit
  character(len=*) :: mesg

  write (unit,'(''Time ['',a19,'']'',a)') f_date_func(), trim(mesg)
  call flush(unit)

end subroutine time_message
  
end module timers
