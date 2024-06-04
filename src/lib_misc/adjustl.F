subroutine adjustl(string)

implicit none

!-pass
character*(*) string
!-local
character(len=1)  :: space
integer, external :: lnblk
integer           :: i, j, k
integer           :: length

length = lnblk(string)

space  = CHAR(32)

if (length.gt.0) then
    i = 1
10  if (string(i:i).eq.space) then
        i = i + 1
        goto 10
    else
        if (i.gt.1) then

            do j=i,length
                k = j - i + 1
                string(k:k) = string(j:j)
            enddo

            do j=length-i+2,length
                 string(j:j) = space
            enddo

        endif
    endif
endif
return
end subroutine
    
