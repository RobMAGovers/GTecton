 character(len=*) function abbrev(cmdname)

 USE MODELDEFINITION

 character(len=*) cmdname
 integer, external :: lnblk

 n = lnblk(cmdname)
 if (cmdname(1:1).eq.'!') then
abbrev = '!'
 else if (n.ge.4 .and. cmdname(1:4).eq.'femi') then
abbrev = 'femi'
 else if (n.ge.2.and. cmdname(1:2).eq.'wf'.and.NTYPE.eq.2) then
abbrev = 'wf'
 else if (n.ge.2 .and. cmdname(1:2).eq.'be') then
abbrev = 'be'
 else if (n.ge.4 .and. cmdname(1:4).eq.'femo') then
abbrev = 'femo'
 else if (n.ge.3 .and. cmdname(1:3).eq.'np ') then
abbrev = 'np'
 else if (n.ge.3 .and. &
   (cmdname(1:3).eq.'gr '.or. cmdname(1:3).eq.'gri') ) then
abbrev = 'gri'
 else if (n.ge.3 .and. cmdname(1:3).eq.'gra') then
abbrev = 'gra'
 else if (n.ge.3 .and. cmdname(1:3).eq.'vol') then
abbrev = 'vol'
 else if (n.ge.4 .and. cmdname(1:4).eq.'stra') then
abbrev = 'stra'
 else if (n.ge.4 .and. cmdname(1:4).eq.'dens') then
abbrev = 'dens'
 else if (n.ge.4 .and. cmdname(1:4).eq.'stre') then
abbrev = 'stre'
 else if (n.ge.2 .and. cmdname(1:2).eq.'ra') then
abbrev = 'ra'
 else if (n.ge.2 .and. cmdname(1:2).eq.'pr') then
abbrev = 'pr'
 else if (n.ge.2 .and. cmdname(1:2).eq.'te') then
abbrev = 'te'
 else if (n.ge.2 .and. cmdname(1:2).eq.'fl') then
abbrev = 'fl'
 else if (n.ge.2 .and. cmdname(1:2).eq.'vi') then
abbrev = 'vi'
 else if (n.ge.2 .and. cmdname(1:2).eq.'di') then
abbrev = 'di'
 else if (n.ge.2 .and. cmdname(1:2).eq.'ma') then
abbrev = 'ma'
 else if (n.ge.5 .and. cmdname(1:5).eq.'estre') then
abbrev = 'estre'
 else if (n.ge.5 .and. cmdname(1:5).eq.'estra') then
abbrev = 'estra'
 else if (n.ge.4 .and. cmdname(1:4).eq.'cons') then
abbrev = 'cons'
 else if (n.ge.3 .and. cmdname(1:3).eq.'era') then
abbrev = 'era'
 else if (n.ge.3 .and. cmdname(1:3).eq.'efl') then
abbrev = 'efl'
 else if (n.ge.3 .and. cmdname(1:3).eq.'edi') then
abbrev = 'edi'
 else if (n.ge.1 .and. cmdname.eq.'?') then
abbrev = '?'
 else if (n.ge.2 .and. cmdname(1:2).eq.'ec') then
abbrev = 'ec'
 else if (n.ge.2 .and. cmdname(1:2).eq.'sh') then
abbrev = 'sh'
 else if (n.ge.2 .and. ( &
   cmdname(1:2).eq.'q ' .or. cmdname(1:2).eq.'qu' .or. &
   cmdname.eq.'q!') ) then
abbrev = 'qu'
 else if (n.ge.1 .and. cmdname(1:1).eq.'#') then
abbrev = '#'
 else if (n.ge.2 .and. cmdname(1:2).eq.'se') then
abbrev = 'se'
 else if (n.ge.2 .and. cmdname(1:2).eq.'mo') then
abbrev = 'mo'
 else if (n.ge.4 .and. cmdname(1:4).eq.'nomo') then
abbrev = 'nomo'
 else if (n.ge.2 .and. cmdname(1:2).eq.'cf') then
abbrev = 'cf'
 else if (n.ge.2 .and. &
 (cmdname(1:2).eq.'e ' .or. cmdname.eq.'e!'))then
abbrev = 'e'
 else if (n.ge.2 .and. cmdname(1:2).eq.'w ') then
abbrev = 'w'
 else if (n.ge.2 .and. cmdname(1:2).eq.'h ') then
abbrev = 'h'
 else if (n.ge.4 .and. cmdname(1:4).eq.'hist') then
abbrev = 'h'
 else if (n.ge.2 .and. cmdname.eq.'we') then
abbrev = 'we'
 else if (n.ge.2 .and. cmdname.eq.'wq') then
abbrev = 'wq'
 else if (n.ge.2 .and. cmdname     .eq.'npl') then
abbrev = 'npl'
 else if (n.ge.2 .and. cmdname(1:2).eq.'r ') then
abbrev = 'r'
 else if (n.ge.2 .and. cmdname(1:2).eq.'sc') then
abbrev = 'sc'
 else if (n.ge.2 .and. cmdname(1:2).eq.'up') then
abbrev = 'up'
 else if (n.ge.4 .and. cmdname(1:4).eq.'elmc') then
abbrev = 'elmc'
 else if (n.ge.4 .and. cmdname(1:4).eq.'surf') then
abbrev = 'surf'
 else if (n.ge.2 .and. cmdname(1:2).eq.'fo') then
abbrev = 'fo'
 else if (n.ge.3 .and. cmdname(1:3).eq.'rot') then
abbrev = 'rot'
 else if (n.ge.4 .and. cmdname(1:4).eq.'conn') then
abbrev = 'conn'
 else if (n.ge.4 .and. cmdname(1:4).eq.'slip') then
abbrev = 'slip'
 else if (n.ge.4 .and. cmdname(1:4).eq.'samp') then
abbrev = 'samp'
 else
if (n.gt.0) then
abbrev = cmdname(1:n)
else
abbrev = ' '
endif
 endif
 return
 end
