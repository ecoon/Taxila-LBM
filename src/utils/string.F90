#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

module String_module
  use petsc
! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private
#include "lbm_definitions.h"

  public :: StringCompare, &
            StringCompareIgnoreCase, &
            StringToUpper, &
            StringToLower, &
            StringReadQuotedWord, &
            StringStartswithAlpha, &
            StringAdjustl, &
            StringNull

  interface StringCompare
    module procedure StringCompare1
    module procedure StringCompare2
  end interface

  interface StringCompareIgnoreCase
    module procedure StringCompareIgnoreCase1
    module procedure StringCompareIgnoreCase2
  end interface

contains

! ************************************************************************** !
!
! StringCompare1: compares two strings
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
PetscBool function StringCompare1(string1,string2,n)

  implicit none

  PetscInt :: i, n
  character(len=n) :: string1, string2
  
  do i=1,n
    if (string1(i:i) /= string2(i:i)) then
      StringCompare1 = PETSC_FALSE
      return
    endif
  enddo

  StringCompare1 = PETSC_TRUE
  return

end function StringCompare1

! ************************************************************************** !
!
! StringCompare2: compares two strings
! author: Glenn Hammond
! date: 10/25/11
!
! ************************************************************************** !
PetscBool function StringCompare2(string1,string2)

  implicit none

  PetscInt :: i, length1, length2
  character(len=*) :: string1, string2
  
  length1 = len_trim(string1)
  length2 = len_trim(string2)
  if (length1 /= length2) then
    StringCompare2 = PETSC_FALSE
    return
  endif

  do i=1,length1
    if (string1(i:i) /= string2(i:i)) then
      StringCompare2 = PETSC_FALSE
      return
    endif
  enddo

  StringCompare2 = PETSC_TRUE
  return

end function StringCompare2

! ************************************************************************** !
!
! StringCompareIgnoreCase1: compares two strings
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
function StringCompareIgnoreCase1(string1,string2,n)

  implicit none

  PetscInt :: i, n
  character(len=n) :: string1, string2
  
  character(len=n) :: upper1, upper2
  PetscBool :: StringCompareIgnoreCase1
  
  upper1 = string1
  upper2 = string2
  
  call StringToUpper(upper1)
  call StringToUpper(upper2)
  
  do i=1,n
    if (upper1(i:i) /= upper2(i:i)) then
      StringCompareIgnoreCase1 = PETSC_FALSE
      return
    endif
  enddo

  StringCompareIgnoreCase1 = PETSC_TRUE
  return

end function StringCompareIgnoreCase1

! ************************************************************************** !
!
! StringCompare: compares two strings
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
function StringCompareIgnoreCase2(string1,string2)

  implicit none

  PetscInt :: i, length1, length2
  character(len=*) :: string1, string2
  
  character(len=MAXSTRINGLENGTH) :: upper1, upper2
  PetscBool :: StringCompareIgnoreCase2
  
  length1 = len_trim(string1)
  length2 = len_trim(string2)
  if (length1 /= length2) then
    StringCompareIgnoreCase2 = PETSC_FALSE
    return
  endif

  upper1 = string1
  upper2 = string2
  
  call StringToUpper(upper1)
  call StringToUpper(upper2)
  
  do i=1,length1
    if (upper1(i:i) /= upper2(i:i)) then
      StringCompareIgnoreCase2 = PETSC_FALSE
      return
    endif
  enddo

  StringCompareIgnoreCase2 = PETSC_TRUE
  return

end function StringCompareIgnoreCase2

! ************************************************************************** !
!
! StringToUpper: converts lowercase characters in a card to uppercase
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine StringToUpper(string)
      
  implicit none

  PetscInt :: i
  character(len=*) :: string

  do i=1,len_trim(string)
    if (string(i:i) >= 'a' .and. string(i:i) <= 'z') then
      string(i:i) = achar(iachar(string(i:i)) - 32)
    endif
  enddo

end subroutine StringToUpper

! ************************************************************************** !
!
! StringToLower: converts uppercase characters in a card to lowercase
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine StringToLower(string)
      
  implicit none

  PetscInt :: i
  character(len=*) :: string

  do i=1,len_trim(string)
    if (string(i:i) >= 'A' .and. string(i:i) <= 'Z') then
      string(i:i) = achar(iachar(string(i:i)) + 32)
    endif
  enddo

end subroutine StringToLower

! ************************************************************************** !
!
! StringReadQuotedWord: reads and removes a name from a string read from the
!                       database.  "'" are used as delimiters.
! author: Glenn Hammond
! date: 11/10/08
!
! ************************************************************************** !
subroutine StringReadQuotedWord(string, name, return_blank_error, ierr)

  implicit none

  PetscInt :: i, begins, ends, realends, length
  PetscBool :: return_blank_error ! Return an error for a blank line
                                ! Therefore, a blank line is not acceptable.
  character(len=*) :: string
  character(len=*) :: name
  PetscBool :: openquotefound
  PetscErrorCode :: ierr

  if (ierr /= 0) return

  openquotefound = PETSC_FALSE
  ! Initialize character string to blank.
  do i=1,len_trim(name)
    name(i:i) = ' '
  enddo

  ierr = 0
  length = len_trim(string)

  ! Remove leading blanks and tabs
  i=1
  do while(string(i:i) == ' ' .or. string(i:i) == achar(9)) 
    i=i+1
  enddo

  if (string(i:i) == "'") then
    openquotefound = PETSC_TRUE
    i=i+1
  endif

  begins=i

  if (openquotefound) then
    do while (string(i:i) /= "'")
      if (i > length) exit
      i=i+1
    enddo
  else
  ! Count # of continuous characters (no blanks, commas, etc. in between)
    do while (string(i:i) /= ' ' .and. string(i:i) /= ',' .and. &
              string(i:i) /= achar(9)) ! 9 = tab
      i=i+1
    enddo
  endif

  realends = i
  ends=i-1

  ! Avoid copying beyond the end of the word (32 characters).
  if (ends-begins > MAXWORDLENGTH - 1) ends = begins + MAXWORDLENGTH - 1

  ! Copy (ends-begins) characters to 'chars'
  name = string(begins:ends)
  ! Remove chars from string
  string = string(realends+1:)

end subroutine StringReadQuotedWord

! ************************************************************************** !
!
! StringStartsWithAlpha: Determines whether a string starts with an alpha char
! author: Glenn Hammond
! date: 10/07/10
!
! ************************************************************************** !
function StringStartsWithAlpha(string)
      
  implicit none

  character(len=*) :: string

  PetscBool :: StringStartsWithAlpha

  string = adjustl(string)

  if ((string(1:1) >= 'a' .and. string(1:1) <= 'z') .or. &
      (string(1:1) >= 'A' .and. string(1:1) <= 'Z')) then
    StringStartsWithAlpha = PETSC_TRUE
  else
    StringStartsWithAlpha = PETSC_FALSE
  endif

end function StringStartsWithAlpha

! ************************************************************************** !
!
! StringAdjustl: Left adjusts a string by removing leading spaces and tabs.
!                This subroutine is needed because the adjustl() Fortran 90 
!                intrinsic will not remove leading tabs.
! author: Richard Tran Mills
! date: 9/21/2010
!
! ************************************************************************** !
subroutine StringAdjustl(string)

  implicit none

  character(len=*) :: string
  
  PetscInt :: i
  PetscInt :: string_length
  character(len=1) :: tab

  ! We have to manually convert any leading tabs into spaces, as the 
  ! adjustl() intrinsic does not eliminate leading tabs.
  tab = achar(9)
  i=1
  string_length = len_trim(string)
  do while((string(i:i) == ' ' .or. string(i:i) == tab) .and. &
           i <= string_length)
    if (string(i:i) == tab) string(i:i) = ' '
    i=i+1
  enddo

  ! adjustl() will do what we want, now that tabs are removed.
  string = adjustl(string) 

end subroutine StringAdjustl

! ************************************************************************** !
!
! StringNull: Returns PETSC_TRUE if a string is blank
! author: Glenn Hammond
! date: 10/07/10
!
! ************************************************************************** !
function StringNull(string)
      
  implicit none

  character(len=*) :: string

  PetscBool :: StringNull
  PetscInt :: length

  length = len_trim(adjustl(string))
  if (length > 0) then
    StringNull = PETSC_FALSE
  else
    StringNull = PETSC_TRUE
  endif

end function StringNull

end module String_module
