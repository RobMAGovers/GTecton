subroutine prerr(mess)
!-pass
character*(*) mess

#ifdef IRIX
call perror(mess)
#endif
#ifdef SOLARIS
call perror(mess)
#endif
#ifdef LINUX
!    call perror(mess)
#endif
    return
    end
