module documentation

implicit none

! This module contains contains several runtime help features.
public  :: plnplt_manual
public  :: plt3d_manual


integer :: plotVariety ! so that the documentation module knows what to display.
                       ! is set by the pltcon subroutine in file [axi/opn/pln]plot.f 
                       ! by values from the constants module.

contains

! Here follow the plnplt and plt3d help function that are called using the question amrtk on the command line


subroutine plnplt_manual(ofwhat)

use constants, only: isAxiplt, isPlnplt, isOpnplt
!use debugmodule, only: plotVariety

character*(*) ofwhat
integer :: helpListWidth, nHelpItems, iHelp


#ifdef EXPERIMENTALTRANSIENTDISCREEP
character(len=14) :: helpList(82)
#else
character(len=14) :: helpList(81)
#endif

helpList = (/"             !",&
             "             #",&
             "             ?",&
             "            be",& 
             "            cf",& 
             "          comp",&
             "          conn",&
             "          cons",&
             "            cp",&
             "          dens",&
             "            di",&
             "            dt",&
             "            du",&
             "             e",&
             "           ebc",&
             "            ec",&
             "           edi",&
             "           efl",&
             "  elementratio",&
             "   elementsize",&
             "          elmc",&
             "           era",&
             "         error",&
             "         estra",&
             "         estre",&
             "           eve",&
             "          femi",&
             "          femo",&
             "            fl",&
             "            fo",&
             "          good",&
             "            gp",&
             "           gra",&
             "          gref",&
             "           gri",&
             "             h",&
             "          heat",&
             "           int",&
             "     invariant",&
             "           lat",&
             "           lon",&
             "            ma",&
             "          mesh",&
             "            mo",&
             "           nbc",&
             "          nomo",&
             "            np",&
             "            pr",&
             "             q",&
             "             r",&
             "            ra",&
             "          read",&
             "        regrid",&
             "          resa",&
             "           rot",&
             "          rotr",&
             "          samp",&
             "          save",&
             "            sc",&
             "            se",&
             "          sedi",&
             "            sh",&
             "          slip",&
             "     slipsides",&
             "          stra",&
             "          stre",&
             "          surf",&
#ifdef EXPERIMENTALTRANSIENTDISCREEP
             "  taylorstress",&
#endif
             "            tc",&
             "            te",&
             "             u",&
             "            ve",&
             "            vi",&
             "           vol",&
             "             w",&
             "          wait",&
             "            we",&
             "            wf",&
             "            wq",&
             "            xl",&
             "            yi",&
             "            yl"/)

if (ofwhat(1:4) .eq. 'list') then
    nHelpItems =  size(helpList)
    helpListWidth=4
    
    write(*,*) 'all available commands:'
    do iHelp=1,nHelpItems-1
        if (mod(iHelp,helpListWidth) .eq. 0) then
            write(*,'(a)'             ) helpList(iHelp)
        else
            write(*,'(a)',advance="no") helpList(iHelp)
        endif
    enddo
    write(*,'(a)'             ) helpList(iHelp)
    write(*,*) ' see ""? [command]"" for help on individual commands'

 else if (ofwhat(1:1) .eq. '!') then
write(*,*) '! command'
write(*,*) '    Issue a shell command'
 else if (ofwhat(1:1) .eq. '#') then
write(*,*) '# comment'
write(*,*) '  records beginning with `#'' will be ignored.'
write(*,*) '  Inline comments start with `\"''. The rest of the record will be ignored.'
 else if (ofwhat(1:1) .eq. '?') then
write(*,*) '? [command]'
write(*,*) '    Print command information. Arguments:'
write(*,*) '        command,        gives verbose info on ""command"".'
 else if (ofwhat(1:2) .eq. 'be') then
write(*,*) 'be[ta] [x|y|z] file'
write(*,*) '    Write extension factor beta as contour data. Flags x, y and z can be'
write(*,*) '    used to select the spatial component of extension (default: x).'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'', `set cst'').'
 else if (ofwhat(1:5) .eq. 'bound') then
write(*,*) 'bound[aries] file'
write(*,*) '    Write line sections describing the partition boundaries to file.'
 else if (ofwhat(1:2) .eq. 'cf') then
write(*,*) 'cf file'
write(*,*) '    Execute commands from command file `file''. (See also: `restart'').'
 else if (ofwhat(1:4) .eq. 'comp') then
write(*,*) 'comp[are] number type InFile OutFile'
write(*,*) '    for timestep ''number'' compare the data ''type'' [totdisp|xdisp|ydisp|i'
write(*,*) '    stress] on nodal points of the 2nd grid read from InFile'
write(*,*) '    (alternative format) with those from the 1st grid (read via femout).'
write(*,*) '    Differences in displacements are written to OutFile (p-type 7).'
 else if (ofwhat(1:4) .eq. 'conn') then
write(*,*) 'conn[ection] file'
write(*,*) '    Writes connectivity file for plotting GMT pscontour data (see ''set gmt'').'
 else if (ofwhat(1:4) .eq. 'cons') then
write(*,*) 'cons[ervation] file'
write(*,*) '    Generates contour data file of mass loss percentages (M - M0) / M0 * 100.'
write(*,*) '    (see also: ""monitor"", ""set"").'
 else if (ofwhat(1:2) .eq. 'cp') then
write(*,*) 'cp FromFile ToFile'
write(*,*) '    copies Fe-disk file FromFile to ToFile, converting binary to ascii'
write(*,*) '    or vice versa.'
 else if (ofwhat(1:4) .eq. 'dens') then
write(*,*) 'dens[ity] file'
write(*,*) '    Write (temperature dependent) densities as contour data.'
write(*,*) '    (See also: `set ascii'', `set deformed'', `monitor'', `set verbose'').'
 else if (ofwhat(1:2) .eq. 'di') then
write(*,*) 'di[splacements] [normal|slip|split|total] file'
write(*,*) '     Writes nodal displacement vector data to file (x, y, angle wrt x-axis, length)'
write(*,*) '     for plotting with GMT "psxy" or "plot".'
write(*,*) '    (See also: `sample'', `set ascii'', `set deformed'', `set npi'', `monitor'','
write(*,*) '    `set arrowhead'', `set verbose'', `set displacement'', `set [no]scaling'').'
 else if (ofwhat(1:2) .eq. 'dt') then
write(*,*) 'dt file [acc|osc|sta] [[no]br[ittle]]'
write(*,*) '    Write time step characteristics to file.'
write(*,*) '        acc[curacy]     relative modal error'
write(*,*) '        osc[illation]   ratio of current time step size to oscillation limit'
write(*,*) '        sta[bility]     ratio of current time step size to stability limit'
write(*,*) '        br[ittle]       also check brittle elements'
write(*,*) '        nobr[ittle]     do not check brittle elements'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'', `set append'', `set cst'').'
 else if (ofwhat(1:2) .eq. 'du') then
write(*,*) 'du[plicate] FromFile ToFile'
write(*,*) '    duplicates FE output file FromFile into ToFile, writing the FE'
write(*,*) '    output in alternative form (X,D,STN).'
! else if ((ofwhat(1:1) .eq. 'e'  .and. len(trim(ofwhat)).eq.1) .or. &
!           ofwhat(1:3) .eq. 'edi') then
!write(*,*) 'e[!] [file]'
!write(*,*) '    edit command file. If ""e!"" use the original file.'
 else if (ofwhat(1:3) .eq. 'ebc') then
write(*,*) 'ebc w|s file'
write(*,*) '    Create files for plotting element mechanical boundary conditions. Flags:'
write(*,*) '    w: Winkler restoring pressure (for psxy)'
write(*,*) '    s: Stresses boundary conditions (for psvelo -Sx)'
write(*,*) '    Flux boundary conditions are not implemented.'
write(*,*) '    In gmt mode no symbol setting and scaling is needed. Set [no]scaling has'
write(*,*) '    influence.'
write(*,*) ' '
write(*,*) '    If set nogmt is on (plotting in P), environment variable PU must be set'
write(*,*) '    to the ""P"" subdirectory putil (e.g. setenv PU /local/P/putil). Flags have'
write(*,*) '    have no influence and can be omitted; all symbols are written to output.'
write(*,*) '    Flux boundary conditions are not implemented. (See also: `nbc'', `set gmt'','
write(*,*) '    `set sysi'', `set deformed'', `set noscaling'').'
 else if (ofwhat(1:2) .eq. 'ec') then
write(*,*) 'ec[ho] [-n] [text]'
write(*,*) '    echoes text. Flag -n suppresses line feed.'
 else if (ofwhat(1:3) .eq. 'edi') then
write(*,*) 'edi[splacements] file [total|normal|slip|split] [e|x|y|z]'
write(*,*) '     Writes nodal displacements magnitudes as contour data to file.'
write(*,*) '     Optional flag determine which displacements are used;'
write(*,*) '         total,         sum of ordinary, slippery node and split node'
write(*,*) '                        displacements'
write(*,*) '         normal,        ordinary nodal displacements'
write(*,*) '         slip,          slippery node displacements'
write(*,*) '         split,         split (i.e. faulted) node displacements'
write(*,*) '         e,             magnitude'
    if (plotVariety.eq.isPlnplt) then
write(*,*) '         x,             x-displacement components'
write(*,*) '         y,             y-displacement components'
else if (plotVariety.eq.isAxiplt) then
write(*,*) '         r,             r-displacement components'
write(*,*) '         z,             z-displacement components'
else if (plotVariety.eq.isAxiplt) then
write(*,*) 'plotting edi in opn not yet tested.'
else
write(*,*) 'plotting variety not recognized.'
write(*,*) 'This should not happen.'
write(*,*) 'Pleaase contact model support.'
endif
write(*,*) ' '
write(*,*) '    When outputting slippery displacements/velocities only,'
write(*,*) '    the sign of the largest displacement/velocity component is added to the magnitude.'
write(*,*) '    This may give alternating sign results for faults with strike azimuths'
write(*,*) '    of about 135 and 315 degrees.'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'').'
 else if (ofwhat(1:3) .eq. 'efl') then
write(*,*) 'efl[ux] file'
write(*,*) '    Write effective heat flux rate as contour data. (See also: `set ascii'','
write(*,*) '    `set deformed'', `set triang'', `monitor'', `set verbose'', `set cst'').'
 else if (ofwhat(1:12) .eq. 'elementratio') then
write(*,*) 'elementratio file'
write(*,*) '    Write quality of elements as contour data.'
write(*,*) '    The quality is defined as four times the area of the inscribed circle'
write(*,*) '    divided by the area of the circumscribed circle. Equilateral triangles have quality 1.'
write(*,*) '    More elongated triangles have a lower quality.'
 else if (ofwhat(1:11) .eq. 'elementsize') then
write(*,*) 'elementsize file'
write(*,*) '    Write area of elements as contour data.'
 else if (ofwhat(1:3) .eq. 'elm') then
write(*,*) 'elmc[enter] file'
write(*,*) '    Write element center coordinates to file. (See also: `set ascii'','
write(*,*) '    `set deformed'', `set elmi'', `monitor'', `set verbose'').'
 else if (ofwhat(1:3) .eq. 'era') then
write(*,*) 'era[te] file [no[log]] [to[tal]|[vi[scous]|pl[astic]] [e|xx|yy|xy|xz|yz]'
write(*,*) '    Write effective strain rate as contour data. There is no default. The last'
write(*,*) '    setting of `rate'' or `erate'' is used to set log. Flags ""total"" and'
write(*,*) '    ""viscous"" can be used to calculate either total or viscous strain rates.'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'', `rate'', `set cst'').'
 else if (ofwhat(1:5) .eq. 'error') then
write(*,*) 'error file'
write(*,*) '    Writes the norm of the a posteriori error vector to a gmt contour file'
write(*,*) '    (See also: `set conn'').' 
 else if (ofwhat(1:5) .eq. 'estra') then

if (plotVariety.eq.isPlnplt) then
write(*,*) 'estra[in] [e|xx|yy|xy|xz|yz] file'
write(*,*) '    Write effective strain as contour data. Arguments are for selecting'
write(*,*) '    strain tensor components.'
write(*,*) '    and xy are stress components, s1 and s3 give largest and smallest eigen-'
write(*,*) '    value. (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'', `set cst'', `set incompressible'').'
else if (plotVariety.eq.isAxiplt) then
write(*,*) 'estra[in] [e|rr|zz|rz|rt|zt|tt] file'
write(*,*) '    Write effective strain as contour data; e=effective stress (default), rr, zz,'
write(*,*) '    rz, rt, zt and tt are strain components.'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'', `set cst'', `set incompressible'').'
else if (plotVariety.eq.isOpnplt) then
write(*,*) ' Help of extre in opnplt not specified. Sorry for that.'
endif


 else if (ofwhat(1:5) .eq. 'estre') then
if (plotVariety.eq.isPlnplt) then
write(*,*) 'estre[ss] [e|xx|yy|xy|s1|s3] file'
write(*,*) '    Write stress-derived contour data; e=effective stress (default), xx, yy,'
write(*,*) '    and xy are stress components, s1 and s3 give largest and smallest eigen-'
write(*,*) '    value. (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'', `set cst'', `set incompressible'').'
else if (plotVariety.eq.isAxiplt) then
write(*,*) 'estre[ss] [e|rr|zz|rz|rt|zt|tt] file'
write(*,*) '    Write stress-derived contour data; e=effective stress (default), rr, zz,'
write(*,*) '    rz, rt, zt and tt are stress components, value. '
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'', `set cst'', `set incompressible'').'
else if (plotVariety.eq.isOpnplt) then
write(*,*) ' Help of extre in opnplt not specified. Sorry for that.'
endif

 else if (ofwhat(1:3) .eq. 'eve') then
write(*,*) 'eve[locity] file [total|normal|slip|split] [e|x|y|z]'
write(*,*) '     Writes nodal velocity magnitudes as contour data to file.'
write(*,*) '     Optional flag determine which velocities are used;'
write(*,*) '         total,         sum of ordinary, slippery node and split node'
write(*,*) '                        velocities'
write(*,*) '         normal,        ordinary nodal velocities'
write(*,*) '         slip,          slippery node velocities'
write(*,*) '         split,         split (i.e. faulted) node velocities'
write(*,*) '         e,             magnitude'
write(*,*) '         x,             x-velocity components'
write(*,*) '         y,             y-velocity components'
write(*,*) '         z,             z-velocity components'
write(*,*) ' '
write(*,*) '    When outputting slippery displacements/velocities only,'
write(*,*) '    the sign of the largest displacement/velocity component is added to the magnitude.'
write(*,*) '    This may give alternating sign results for faults with strike azimuths'
write(*,*) '    of about 135 and 315 degrees.'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'').'
 else if (ofwhat(1:4) .eq. 'femi') then
write(*,*) 'femi[n] [MechFile [TempFile]]'
write(*,*) '    Read finite element (TECTON) input file(s). ""MechFile"" (typically named'
write(*,*) '    ""TECIN.DAT"") contains mesh geometry, material properties and boundary'
write(*,*) '    conditions for the mechanical problem. If neccesary, ""TempFile"" (typically'
write(*,*) '    named ""TEMPIN.DAT"") is read. ""TempFile"" contains material properties and'
write(*,*) '    boundary conditions relevant for the thermal problem. (See also: `femo'').'
 else if (ofwhat(1:4) .eq. 'femo') then
write(*,*) 'femo[ut] [file] [bin[ary]|asc[ii]]'
write(*,*) '    Read finite element output files. ""file"" is the generic finite element output'
write(*,*) '    file which should be similar to the output file specified for pln being of the'
write(*,*) '    form ""[dirname/][basename]"". If nothing is specified the default ""FEDSK.DAT"" is'
write(*,*) '    used. Normally the input format (binary or ascii) is determined automatically.'
write(*,*) '    The file format can be forced by giving ""binary"" or ""ascii"" arguments.'
write(*,*) '    (See also: `femi'').'
 else if (ofwhat(1:2) .eq. 'fl') then
write(*,*) 'fl[ux] file'
write(*,*) '    Write element heat fluxes to file as arrows. (See also: `set ascii'','
write(*,*) '    `set deformed'', `set elmi'', `monitor'', `set arrowhead'', `set verbose'','
write(*,*) '    `set flux'', `set [no]scaling'').'
 else if (ofwhat(1:2) .eq. 'fo') then
write(*,*) 'fo[rce] [file]'
write(*,*) '    Writes information on force level in monitor elements.'
 else if (ofwhat(1:4) .eq. 'good') then
write(*,*) 'good[lement] file'
 else if (ofwhat(1:2) .eq. 'gp') then
write(*,*) 'gp lon lat azimuth'
write(*,*) '        Sets equator for geographic projection. (See also: ""set project"",'
write(*,*) '        ""lon"", ""lat"", ""xl"", ""yl"", ""pr"").'
 else if (ofwhat(1:3) .eq. 'gra') then
write(*,*) 'gra[vity] [xstep [xmin xmax]] file'
write(*,*) '    Calculate gravity anomalies along surface. Should be preceded by definition'
write(*,*) '    of reference gravity g0. The raw anomaly is defined by'
write(*,*) ' '
write(*,*) '                dg(x,y) = g(x,y) - g0'
write(*,*) ' '
write(*,*) '    where x and y are model upper surface coordinates. In this definition, a mass deficit'
write(*,*) '    produces a positive anomaly, and a mass excess is related with a negative anomaly.'
write(*,*) ' '
write(*,*) '    Subsequent gravity reductions are controlled by the setting of the'
write(*,*) '    gravity mode:'
write(*,*) '        If ""raw"", no reduction occurs.'
write(*,*) '        If ""free-air"",'
write(*,*) ' '
write(*,*) '                dg(x,y) = g(x,y) - g0 + dgh'
write(*,*) ' '
write(*,*) '        where, dgh is the reduction to reference level Yref.'
write(*,*) ' '
write(*,*) '        If ""bouguer"",'
write(*,*) ' '
write(*,*) '                dg(x,y) = g''(x,y) - g0'
write(*,*) ' '
write(*,*) '        Bouguer anomalies are calculated by integration of gravity contributions from'
write(*,*) '        below reference level Yref, i.e. not by reducing the raw anomaly! Although'
write(*,*) '        this is more accurate for calculating the real Bouguer anomaly, results may appear'
write(*,*) '        inconsistent with free-air and raw anomalies.'
write(*,*) ' '
write(*,*) '    (See also: `set expand'', `gref'', `set bouguer'', `set raw'', `set free-air'','
write(*,*) '    `set density'') .'
 else if (ofwhat(1:4) .eq. 'gref') then
write(*,*) 'gref x'
write(*,*) '    Calculate (raw)  reference gravity signal g0 at X=x, Y=Yref.'
write(*,*) '    Yref is set to 0 initially, and may be reset using either the'
write(*,*) '    `set bouguer Yref'' or `set free-air Yref'' commands.'
write(*,*) '    (See also: ""grav"", ""set expand"", ""set density"", ""set bouguer"",'
write(*,*) '    ""set free-air"")'
 else if (ofwhat(1:2) .eq. 'gr') then
write(*,*) 'gr[id] file [[no]nu[mbers]]'
write(*,*) '    Generates "P" data- or command-file of the finite element mesh.'
write(*,*) '    Output will be written to `file''. In case argument `numbers'' is given'
write(*,*) '    a command file with `pst'' is written. There is no default. Last setting'
write(*,*) '    is used. (See also: `set deformed'', `set ascii'', `set append'', `monitor'','
write(*,*) '    `np'', `set shrink'') '
 else if (ofwhat(1:1) .eq. 'h') then
write(*,*) 'h[istory] [number]'
write(*,*) '    prints history of last commands.'
 else if (ofwhat(1:4) .eq. 'heat') then
write(*,*) 'heat file'
write(*,*) '    Write shear heating rate per unit volume as contour data.'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'', `set cst'').'
 else if (ofwhat(1:3) .eq. 'int') then
write(*,*) 'int[erpolate] number type InFile OutFile'
write(*,*) '    for timestep ''number'' interpolate nodal points of the 2nd grid read'
write(*,*) '    from InFile (alternative format) on 1st grid (read via femout).'
write(*,*) '    Interpolated displacements are written to OutFile (p-type 7)'
write(*,*) '    Datatypes : [tot|x|y]'
 else if (ofwhat(1:9) .eq. 'invariant') then
write(*,*) 'invariant []1/2/3] outFile'
write(*,*) '    writes first, second or third invariant of the stress tensor to outFile.'
 else if (ofwhat(1:3) .eq. 'lat') then
write(*,*) 'lat latmin dlat latmax'
write(*,*) '        Sets latitude window for geographic projection. (See also:'
write(*,*) '        ""set project"", ""gp"", ""lon"", ""xl"", ""yl"", ""pr"").'
 else if (ofwhat(1:3) .eq. 'lon') then
write(*,*) 'lon lonmin dlon lonmax'
write(*,*) '        Sets longitude window for geographic projection. (See also:'
write(*,*) '        ""set project"", ""gp"", ""lat"", ""xl"", ""yl"", ""pr"").'
 else if (ofwhat(1:2) .eq. 'ma') then
write(*,*) 'ma[terial] [flag] file'
write(*,*) '    Writes material property data to file'
write(*,*) '    ""flag"" should be one of:'
write(*,*) '        nu[mber],       material number'
write(*,*) '        E,              Young''s modulus'
write(*,*) '        v,              Poisson''s ratio'
write(*,*) '        vi[scosity]     viscosity'
write(*,*) '        po[wer],        viscosity stress power'
write(*,*) '        de[nsity],      density'
write(*,*) '        th[ick],        element thickness'
write(*,*) '        Qpl,            powerlaw activation energy'
write(*,*) '        Apl,            powerlaw pre-exponential constant'
write(*,*) '        npl,            stress power for powerlaw creep'
write(*,*) '        Vpl,            powerlaw activation volume'
write(*,*) '        pl[astype]      plasticity type'
write(*,*) '        fr[ictangle]    plasticity friction angle'
write(*,*) '        un[iyld]        (equivalent) uniaxial yield stress'
write(*,*) '        fl[uidity]      fluidity'
write(*,*) '        yld[pwr]        plastic yield power'
write(*,*) '        nflow           plastic yielding law'
write(*,*) '        ha[rding]       plastic hardening'
write(*,*) '        he[at],         heat production'
write(*,*) '        Cp,             specific heat'
write(*,*) '        al[pha],        thermal expansion'
write(*,*) '        Kx,             X-conductivity'
write(*,*) '        Ky,             Y-conductivity'
write(*,*) '    (See also: `set deformed'', `monitor'', `set verbose'', `viscosity'').'
 else if (ofwhat(1:4) .eq. 'mesh') then
write(*,*) 'mesh [informationType] [file]'
write(*,*) '    The mesh command displays information about the mesh and writes this'
write(*,*) '    as contour data.'
write(*,*) '    the information type can be:'
write(*,*) '    "area",        Writing the area of the triangle'
write(*,*) '    "arearatio",   Writing the element quality, defined as four times '
write(*,*) '                   the area of the '
write(*,*) '                   inscribed circle divided by the area of the circumscribed '
write(*,*) '                   circle. Equilateral triangles have quality 1.'
write(*,*) '                   More elongated triangles have a lower quality.'
write(*,*) '    "sideratio",   Writing the element quality, defined as the length of the'
write(*,*) '                   shortest side divided by the length of the longest side.'
 else if (ofwhat(1:2) .eq. 'mo') then
write(*,*) 'mo[nitor] flag number[s]'
write(*,*) '    monitor output at a specific node or element number. There is no'
write(*,*) '    default. Turned off by `nomonitor nodes'' or `nomonitor elements''.'
write(*,*) '        This command affects the syntax of commands that normally require a'
write(*,*) '    file name argument. If, in these commands, the file name argument is'
write(*,*) '    omitted, information is written to standard output. If the file name'
write(*,*) '    argument is given, monitor data are written (ASCII format) to this'
write(*,*) '    file. If `noverbose'' has been set, monitor information is very compact.'
write(*,*) '        ""flag"" should be one of:'
write(*,*) '            n[odes],            NUMBERS are interpreted as nodal point numberss.'
write(*,*) '            e[lements],         NUMBERS are interpreted as element numners.'
write(*,*) '    Ranges of numbers can be given as ""istart,istep,iend"" or ""istart,iend""'
write(*,*) '    A valid statement is therefore:'
write(*,*) '        monitor e 3,5,28'
write(*,*) '    (See also: `set verbose'').'
! else if (ofwhat .eq. 'nbc' .and. mattsym.eq.0) then
! write(*,*) 'nbc [a|d|v|f|i|r] file'
! write(*,*) '    Write data to file that allow plotting of nodal boundary condition (bc) '
! write(*,*) '    symbols.'
! write(*,*) ' '
! write(*,*) '    If the output is for plotting with GMT (see `set [no]gmt''), optional'
! write(*,*) '    command flags have the the following meaning:'
! write(*,*) '        a: all bc types'
! write(*,*) '        d: displacement bc'
! write(*,*) '        v: velocity bc'
! write(*,*) '        f: force bc'
! write(*,*) '        i: initial displacement bc'
! write(*,*) '        r: strain rate bc'
! write(*,*) '    Output suitable for plotting with ""psxy -Sv""'
! write(*,*) '    Temperature bcs have not been implemented.'
! write(*,*) '    Symbol size definitions (`set sysi'') do not influence the scaling.'
! write(*,*) ' '
! write(*,*) '    If the output is for plotting with P (see `set [no]gmt''), the'
! write(*,*) '    environment variable PU needs to be set to the P subdirectory putil'
! write(*,*) '    (e.g. setenv PU /local/P/putil). The constant strain rate bc and'
! write(*,*) '    temperature bc have not been implemented. (See also: `ebc'', `set ascii'','
! write(*,*) '    `set append'', `set mattsym'', `set deformed'', `set [no]scaling'').'
 else if (ofwhat(1:3) .eq. 'nbc') then
write(*,*) 'nbc [a|d|v|f|i|r] file'
write(*,*) '    Write data to file that allow plotting of nodal boundary condition (bc) '
write(*,*) '    symbols with GMT'
write(*,*) ' '
write(*,*) '    Flags have the the following meaning:'
write(*,*) '        a: all bc types'
write(*,*) '        d: displacement bc'
write(*,*) '        v: velocity bc'
write(*,*) '        f: force bc'
write(*,*) '        i: initial displacement bc'
write(*,*) '        r: strain rate bc'
write(*,*) '    Output suitable for plotting with "psxy:"'
write(*,*) '            gmt psxy bc_file -J -R -Skbcs -K -O >> bcs.ps'
write(*,*) ' '
write(*,*) '    Temperature bcs have not yet been implemented.'
write(*,*) ' '
write(*,*) '    (See also: `ebc'',`set [no]scaling'',`set mattsym'').'
 else if (ofwhat(1:4) .eq. 'nomo') then
write(*,*) 'nomo[nitor] n[ode]|e[lement]|r[ow]|c[olumn]'
write(*,*) '    Turn monitoring of nodes, elements, rows or columns off.'
 else if (ofwhat(1:2) .eq. 'np') then
write(*,*) 'np file [[no]numbers] [tecton] [[av]erage|[to]tal] / [raw]'
write(*,*) '    Write nodal point coordinates as a ""P"" type 2 data to `file'' in case'
write(*,*) '    of nonumbers. Write ""P"" command file in case numbers has been set.'
write(*,*) '    ""average"" or ""total"" can be used in case that slippery or split nodes'
write(*,*) '    need to be plotted; ""average"" will result in nodal point mean coordinate,'
write(*,*) '    ""total"" will add both slippery and split to the coordinates. Note: correct'
write(*,*) '    ordering of nodal point coordinates is not guaranteed when ""total"" is on.'
write(*,*) '    For all settings there is no default. Last setting is used.'
write(*,*) '      With the ""tecton"" flag, nodal points are written in TECTON nodal point'
write(*,*) '    coordinate input format, e.g. for restart purposes. This flag does not work'
write(*,*) '    when slippery nodes are present and is always false by default.'
write(*,*) '    (See also: `grid'').'
write(*,*) ' '
write(*,*) '    Using only the filename and the ""raw"" flag, a file is written'
write(*,*) '    containing only x y pairs for the coordinates of the nodal points.'
 else if (ofwhat(1:4) .eq. 'part') then
write(*,*) 'part[ition] file'
write(*,*) '    Write a file that contains the partition of each element, so that the'
write(*,*) '    partitioning may be plotted, using GMT'
 else if (ofwhat(1:2) .eq. 'pr') then
write(*,*) 'pr[essure] file'
write(*,*) '    Write pressures as contour data. (See also: `set ascii'','
write(*,*) '    `set deformed'', `set triang'', `monitor'', `set verbose'', `set cst'').'
 else if (ofwhat(1:1) .eq. 'q') then
write(*,*) 'q[uit]'
write(*,*) '    stops execution of program'
 else if (ofwhat(1:1) .eq. 'r' .and. len(trim(ofwhat)).eq.1) then
write(*,*) 'r'
write(*,*) '    restart program. If ""r!"" restart from original command file. Most variable'
write(*,*) '    settings are remembered. (See also: `cf'').'
 else if (ofwhat(1:2) .eq. 'ra') then
write(*,*) 'ra[te] file [no[log]] [to[tal]|vi[scous]|pl[astic]]'
write(*,*) '    Write file for plotting element strain rates as arrows in principal stress'
write(*,*) '    directions using psvelo -Sx. There is no default. The last setting of `rate'' '
write(*,*) '    or `erate'' is used to set logarithmic scaling or not. Flags ""total"" and'
write(*,*) '    ""viscous"" can be used to calculate either total or viscous strain rates.'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set elmi'', `monitor'', `set'' '
write(*,*) '    arrowhead'', `set verbose'', `era'', `set rate'', `set [no]scaling''). '
 else if (ofwhat(1:4) .eq. 'read') then
write(*,*) 'read file'
write(*,*) '    Read nodal points file. Used for resampling nodal point data.'
write(*,*) '    (see also: "resample").'
 else if (ofwhat(1:6) .eq. 'regrid') then
write(*,*) 'regrid material_number[s]'
write(*,*) '    Dynamic Lagrangian Remeshing facility, consisting of the following steps;'
write(*,*) '    1) Reads ""RESTART.DAT"" file.'
write(*,*) '    2) Creates a new Delauney tessellation of nodal points in a contiguous'
write(*,*) '       region with specified material numbers (using system call'
write(*,*) '       ""triangle -pQ /tmp/triangle.poly"").'
write(*,*) '    3) Maps element quantities (stress, strain and slippery displacements) onto'
write(*,*) '       the new grid.'
write(*,*) '    4) Overwrite ""RESTART.DAT"" '
write(*,*) '    5) Write new tessellation file ""tecin.dat.elm_reg""'
write(*,*) '    6) In case of slippery nodes, write ""tecin.dat.sns_reg""'
 else if (ofwhat(1:4) .eq. 'resa') then
write(*,*) 'resa[mple] file'
write(*,*) '    Writes resampled nodal temperature data in tecton format'
write(*,*) '    to file. (see also: ""read"")'
 else if (ofwhat(1:4) .eq. 'rot ' .or. ofwhat(1:4).eq.'rota') then
write(*,*) 'rot[ation] file'
write(*,*) '    Writes antisymmetric part of the displacement gradient tensor to file.'
write(*,*) '    More specifically, in two dimensions, the only nonzero component of this'
write(*,*) '    infinitesimal rotation tensor is 0.5*( dUy/dx - dUx/dy ).'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'', `set cst'').'
 else if (ofwhat(1:4) .eq. 'rotr') then
write(*,*) 'rotr[ate] file'
write(*,*) '    Writes antisymmetric part of the velocity gradient tensor to file.'
write(*,*) '    More specifically, in two dimensions, the only nonzero component of this'
write(*,*) '    infinitesimal rotation rate tensor is 0.5*( dVy/dx - dVx/dy ).'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'', `set cst'').'
 else if (ofwhat(1:4) .eq. 'samp') then
write(*,*) 'samp[le] d|v|T|s|e|c|w|r|f coordinate[s] [filename]'
write(*,*) '   compute quantities at coordinate[s]'
write(*,*) '   d = displacement'
write(*,*) '   v = velocity'
write(*,*) '   T = temperature'
write(*,*) '   s = stress'
write(*,*) '   e = strain'
write(*,*) '   c = viscosity'
write(*,*) '   w = rotation'
write(*,*) '   r = strain rate'
write(*,*) '   f = heat flux'
write(*,*) '   '
write(*,*) '   Data can be sampled in the following ways:'
write(*,*) '   '
write(*,*) '   1: sample d|v|T|s|e|c|w|r|f <xrange> <yrange>'
write(*,*) '   2: sample d|v|T|s|e|c|w|r|f <xrange> <yrange> filename'
write(*,*) '   3: sample d|v|T|s|e|c|w|r|f <inputfile> filename'
write(*,*) '   4: sample d|v|T|s|e|c|w|r|f line <xstart,ystart> <xend,yend> filename'
write(*,*) '   5: sample d|v|T|s|e|c|w|r|f line <xstart,ystart> <xend,yend> [npoints] filename'
write(*,*) '   '
write(*,*) '   A coordinate range can defined as:'
write(*,*) '   '  
write(*,*) '   x'
write(*,*) '   xmin, xmax'
write(*,*) '   xmin, xmax, xincrement'
write(*,*) '   '  
write(*,*) '   and similar for y'
write(*,*) '   '  
write(*,*) '   for the case of: '
write(*,*) '   xmin, xmax'
write(*,*) '   the increment is computed automatically, as follows, in pseudocode: '
write(*,*) '   in order to get round numbers in the magnitude.'
write(*,*) '   '
write(*,*) '   |   range = xman - xmin'
write(*,*) '   |   pow   = 10log(range)'
write(*,*) '   |   ipow  = int(pow)'
write(*,*) '   |   diff  = pow - ipow '
write(*,*) '   |'
write(*,*) '   |        if diff <= 0.176 then    pre =  2'
write(*,*) '   |   else if diff <= 0.398 then    pre =  5'
write(*,*) '   |   else if diff <= 0.778 then    pre = 10'
write(*,*) '   |   else                          pre = 20'
write(*,*) '   |'
write(*,*) '   |   if      range >= 1    then nSteps = pre*10^(ipow-1)'
write(*,*) '   |   else                       nSteps = pre*10^(ipow-2)'
write(*,*) '   |'
write(*,*) '   |   xinc  = (xman - xmin) / (nSteps-1)'
write(*,*) '   '
write(*,*) '   Sample is specifically used to interpolate (based on shape functions) nodal'
write(*,*) '   quantities. Use ""monitor np"" in combination with ""disp"", ""vel"" or ""temp"" to'
write(*,*) '   get nodal values. For element quantities, functionality of sample partly'
write(*,*) '   overlaps with that of ""monitor e"" in combination with ""stress"", ""strain"",'
write(*,*) '   ""rate"" or ""flux"".'
write(*,*) '   '
write(*,*) '   Note that the name of input file with a list of coordinate pairs should not be ''line''.'
write(*,*) '   That will confuse poor little plt3d beyond reasaon.'
write(*,*) '   '
 else if (ofwhat(1:4) .eq. 'save') then
write(*,*) 'save'
write(*,*) '    store current erosion surface.'
 else if (ofwhat(1:2) .eq. 'sc') then
write(*,*) 'sc[alebar] flag x y file [h[orizontal]|v[ertical]]'
write(*,*) '    Writes coordinates of an arrow centered at (x,y) to display the current'
write(*,*) '    scale bar size. Flag is one of'
write(*,*) '        di[splacement]: displacement scale bar'
write(*,*) '        ve[locity]:     velocity scale bar'
write(*,*) '        stre[ss]:       stress scale bar'
write(*,*) '        stra[in]:       strain scale bar'
write(*,*) '        ra[te]:         strain rate scale bar'
write(*,*) '        fl[ux]:         heat flux scale bar.'
write(*,*) '    Last argument controls whether the arrow will be horizontal or vertical.'
write(*,*) '    The ""file""-argument is used as a seed to generate two files:'
write(*,*) '        file.dat:       arrow data'
write(*,*) '        file.cf:        ''pt'' command file'
write(*,*) '    (See also: `set arrowhead'', `set append'', `set ascii'').'
 else if ((ofwhat(1:3) .eq. 'se ').or.(ofwhat(1:3) .eq. 'set')) then
write(*,*) 'se[t] parameter [value[s]]'
write(*,*) ' '
! write(*,*) '    se[t] ar[rowhead] size angle [no[co[nt]]]'
! write(*,*) '        Define size and top angle of arrows.'
! write(*,*) ' '
write(*,*) '    se[t] [no]ap[pend]'
write(*,*) '        Append data to file.'
write(*,*) ' '
write(*,*) '    se[t] [no]as[cii]'
write(*,*) '        Output file format. GMT output is always ASCII'
write(*,*) ' '
! write(*,*) '    se[t] bcsy[mbol] type number'
! write(*,*) '        For P; set symbol number of boundary condition of type `type''.'
! write(*,*) '        `type'' can be one of;'
! write(*,*) '        fi[xed]:        for displacements fixed in more than one direction'
! write(*,*) '        sl[ide]:        for displacements fixed in one direction'
! write(*,*) '        fo[rce]:        for nodal force boundary conditions'
! write(*,*) '        ve[locity]:     for nodal velocity boundary conditions'
! write(*,*) '        st[ress]:       for element stress boundary conditions'
! write(*,*) '        wi[nkler]:      for nodal Winkler boundary conditions'
! write(*,*) ' '
write(*,*) '    se[t] [no]bi[nary]'
write(*,*) '        Output file format. GMT output is always ASCII.'
write(*,*) ' '
write(*,*) '    se[t] bou[guer] [Yref]'
write(*,*) '        Sets gravity calculation mode to ""bouguer"". In this mode, the'
write(*,*) '        gravity attraction of masses below the horizon Y=Yref are calculated.'
write(*,*) '        Yref is set to 0 initially, and may be reset by both ""set bouguer""'
write(*,*) '        and ""set free-air"" commands. (See also ""set raw"", set ""free-air"",'
write(*,*) '        ""gref"", ""grav"", ""set expand"", ""set density"").'
write(*,*) ' '
write(*,*) '    se[t] [no]cst'
write(*,*) '        Sets or unsets Constant STress mode for interpolation. (See also ""set'
write(*,*) '        triangulate"").'
write(*,*) ' '
write(*,*) '    se[t] [no]de[formed] [amplification]'
write(*,*) '        Original or deformed mesh output. Amplification factor can be used'
write(*,*) '        to exaggerate deformation.'
write(*,*) ' '
write(*,*) '    se[t] den[sity] number'
write(*,*) '        Sets density of elements which are y-expanded to a depth of 300 km.'
write(*,*) '        This option is used for calculating gravity anomalies. If the density'
write(*,*) '        is set to a value less or equal zero, current elements densities'
write(*,*) '        will be used (See also ""set raw"", set ""free-air"", ""set bouguer"",'
write(*,*) '        ""gref"", ""grav"", ""set expand"").'
write(*,*) ' '
! write(*,*) '    se[t] di[splacement] value1 [unit1 = value2 unit2]'
! write(*,*) '        Functionality of this command depends on whether output is for GMT'
! write(*,*) '        or for P. For GMT, it defines a multiplication factor on displacement'
! write(*,*) '        data. For P, this command defines the scaling of displacements'
! write(*,*) '        to grid units. ""unit""-fields are ignored on input, but are'
! write(*,*) '        convenient to keep track of the meaning of numbers. Example:'
! write(*,*) '                set displacement 10 m = 5000 m'
! write(*,*) '        The first argument is also used to define the size of a scale bar.'
! write(*,*) '        (see also: `set gmt'', `set [no]scaling'') '
! write(*,*) ' '
write(*,*) '    se[t] [no]ec[ho] [number]'
write(*,*) '        Set echoing mode. `noecho'' is identical to `echo 0''.'
write(*,*) '        echo 1:   verbose mode.'
write(*,*) '        echo 2:   input files reading.'
write(*,*) '        echo 3,4: program calculations'
write(*,*) '        echo 5:   verbose program calculations'
write(*,*) ' '
write(*,*) '    se[t] elmi[nterval] [base] number'
write(*,*) '        Sets element base and increment for plotting of arrows.'
write(*,*) ' '
write(*,*) '    se[t] ex[pand] x|y iside numbers'
write(*,*) '        Sets element numbers and element side numbers which have to be'
write(*,*) '        expanded horizontally or vertically during gravity calculations.'
write(*,*) '        Expansion is done by 10 times maximum x- or y-values.'
write(*,*) '        (See also ""gref"", ""grav"", ""set raw"", ""set bouguer"", ""set free-air"").'
write(*,*) ' '
! write(*,*) '    se[t] fl[ux] value1 [unit1 = value2 unit2]'
! write(*,*) '        Functionality of this command depends on whether output is for GMT'
! write(*,*) '        or for P. For GMT, it defines a multiplication factor on heat flux'
! write(*,*) '        data. For P, this command defines the scaling of fluxes'
! write(*,*) '        to grid units. ""unit""-fields are ignored on input, but are'
! write(*,*) '        convenient to keep track of the meaning of numbers. Example:'
! write(*,*) '                set flux 60e-3 W/m^2 = 5000 m'
! write(*,*) '        The first argument is also used to define the size of a scale bar.'
! write(*,*) '        (see also: `set gmt'', `set [no]scaling'') '
! write(*,*) ' '
write(*,*) '    se[t] fre[e-air] [Yref]'
write(*,*) '        Sets gravity calculation mode to ""free-air"". In this mode,'
write(*,*) '        the gravity attraction of masses below the surface are calculated.'
write(*,*) '        Subsequently, the signal is reduced to the horizon Y=Yref.'
write(*,*) '        Yref is set to 0 initially, and may be reset by both ""set bouguer"" '
write(*,*) '        and ""set free-air"" commands. (See also ""set bouguer"", ""set raw"",'
write(*,*) '        ""gref"", ""grav"", ""set expand"", ""set density"").'
write(*,*) ' '
! write(*,*) '    se[t] [no]gmt'
! write(*,*) '        Write plot output in GMT format (this is the default). Alternatively,'
! write(*,*) '        output will be in P format. When binary output is possible and set'
! write(*,*) '        ("set bin"), 4-byte reals and integers are written.'
! write(*,*) ' '
write(*,*) '    se[t] in[compressible]'
write(*,*) '        Switch for calculating effective stress based on incompressible or'
write(*,*) '        compressible assumption.'
write(*,*) '        Default = False; compressible.'
write(*,*) '        See also: ''set comp'''
write(*,*) ' '
write(*,*) '    se[t] it[ime] number'
write(*,*) '        Set output time step number. Increment ""++"" and decrement ""--""'
write(*,*) '        operators are also allowed.'
write(*,*) ' '
write(*,*) '    se[t] [no]mattsym [c|g]'
write(*,*) '        Set GMT custom boundary condition symbols on [off]. These are designed to be'
write(*,*) '        used with GMT 5.2.1. After turning mattsym option on, a GMT custom symbol'
write(*,*) '        macro file called "bcs.def" will be created when boundary conditions are'
write(*,*) '        requested from plnplt. Plot these custom boundary conditions with psxy:'
write(*,*) '            gmt psxy bc_file -J -R -Skbcs -K -O >> bcs.ps'
write(*,*) '        Optional "c" and "g" arguments are short for "cartesian" and "geographic".'
write(*,*) '        Cartesian [default] for x-y projection, and g for geographic projections'
write(*,*) ' '
write(*,*) '    se[t] npi[nterval] [base] number'
write(*,*) '        Sets nodal point base and increment for plotting of arrows.'
write(*,*) ' '
write(*,*) '    se[t] opn [flag [NUMZ [scale [axnu]]]]'
write(*,*) '        Option for ""opnplt"". Sets ""out-of-plane"" mode. These settings define'
write(*,*) '        the way opn fe-input files are converted to f3d files. ""flag""'
write(*,*) '        determines whether grid dimensions (""x"") or displacements (""d"") are'
write(*,*) '        used to specify out of plane dimensions. NUMNZ specifies the number of'
write(*,*) '        additional nodal point planes which is being generated. ""scale"" is a'
write(*,*) '        linear scaling parameter. ""axnu"" specifies the fraction of the derived'
write(*,*) '        dimensions that is being added to the maximum and subtracted from the'
write(*,*) '        minimum dimension values. (see also: wfemi) '
write(*,*) '        Example: set opn d 4 1.0 0.05'
write(*,*) ' '
write(*,*) '    se[t] [no]pr[oject]'
write(*,*) '        Switches geographic projection off (""nopr"") or on (""pr""). This allows'
write(*,*) '        eigenvectors to be computed in a cm-frame rather than in a longitude-'
write(*,*) '        latitude frame. Consequentely, eigenvectors will be plotted in'
write(*,*) '        perpendicular. Downside is that the equator needs to be specified with'
write(*,*) '        the ""gp"" command, longitude and latitude ranges with ""lon"" and ""lat""'
write(*,*) '        commands, plot physical dimensions with ""xl"" and ""yl"" commands.'
write(*,*) '        In the future, projection type also will need to be specified using'
write(*,*) '        the ""pr"" command. Currently, only the cartesian projection of ""P"" is'
write(*,*) '        implemented.'
write(*,*) ' '
! write(*,*) '    se[t] ra[te] value1 [unit1 = value2 unit2]'
! write(*,*) '        Functionality of this command depends on whether output is for GMT'
! write(*,*) '        or for P. For GMT, it defines a multiplication factor on strain rate'
! write(*,*) '        data. For P, this command defines the scaling of strain rates'
! write(*,*) '        to grid units. ""unit""-fields are ignored on input, but are'
! write(*,*) '        convenient to keep track of the meaning of numbers. Example:'
! write(*,*) '                set rate 1e-14 /s = 5000 m'
! write(*,*) '        The first argument is also used to define the size of a scale bar.'
! write(*,*) '        (see also: `set gmt'', `set [no]scaling'') '
! write(*,*) ' '
write(*,*) '    se[t] [no]sc[aling]'
write(*,*) '        Turns symbol automatic scaling on or off.'
write(*,*) ' '
write(*,*) '    se[t] raw'
write(*,*) '        Sets gravity calculation mode to ""raw"". Gravity calculations, performed'
write(*,*) '        in this mode will result in the gravity anomaly as if measured along'
write(*,*) '        the current model surface. (See also ""set bouguer"", ""set free-air"",'
write(*,*) '        ""grav"", ""gref"", ""set expand"", ""set density"").'
write(*,*) ' '
write(*,*) '    se[t] sh[rink] [factor]'
write(*,*) '        Shrinks elements by a factor. Useful for testing nodal point'
write(*,*) '        connectivity. (See also ""grid"").'
write(*,*) ' '
! write(*,*) '    se[t] stra[in] value [unit1 = value2 unit2]'
! write(*,*) '        Functionality of this command depends on whether output is for GMT'
! write(*,*) '        or for P. For GMT, it defines a multiplication factor on strain'
! write(*,*) '        data. For P, this command defines the scaling of strain'
! write(*,*) '        to grid units. ""unit""-fields are ignored on input, but are'
! write(*,*) '        convenient to keep track of the meaning of numbers. Example:'
! write(*,*) '                set strain 0.01 unit = 5000 m'
! write(*,*) '        The first argument is also used to define the size of a scale bar.'
! write(*,*) '        (see also: `set gmt'', `set [no]scaling'') '
! write(*,*) ' '
! write(*,*) '    se[t] stre[ss] value [comment = unit comment]'
! write(*,*) '        Functionality of this command depends on whether output is for GMT'
! write(*,*) '        or for P. For GMT, it defines a multiplication factor on stress'
! write(*,*) '        data. For P, this command defines the scaling of stresses'
! write(*,*) '        to grid units. ""comment""-fields are ignored on input, but are'
! write(*,*) '        convenient to keep track of the meaning of numbers. Example:'
! write(*,*) '                set stress 2.E8 Pa = 5.E3 m'
! write(*,*) '        The first argument is also used to define the size of a scale bar.'
! write(*,*) '        (see also: `set gmt'', `set [no]scaling'') '
! write(*,*) ' '
! write(*,*) '    se[t] sy[mbolsize] size'
! write(*,*) '        Set boundary condition symbol size in plot units.'
! write(*,*) ' '
write(*,*) '    se[t] [no]sign'
write(*,*) '        Add sign information to tensor magnitudes (effective stress, strain'
write(*,*) '        and so on). This only works in two dimensions.'
write(*,*) ' '
write(*,*) '    se[t] noso[elim] or dosoelim'
write(*,*) ' '
write(*,*) '    se[t] tem[perature] a[bsolute]|r[elative]|i[ncrement] [itime]'
write(*,*) '        Set temperatures mode.'
write(*,*) '            `absolute''  =current temperatures'
write(*,*) '            `relative''  =temperature differences relative to some other'
write(*,*) '                         time, specified by itime (default: itime=0) '
write(*,*) '            `increment'' =temperature differences relative to previous output'
write(*,*) ' '
write(*,*) '    se[t] [no]test'
write(*,*) '        If set, additional tests on the GTECTON input data will be performed'
write(*,*) '        Default=notest'
write(*,*) ' '
write(*,*) '    se[t] [no]tr[iangulate]'
write(*,*) '        Triangulation mode. Usefull when contouring on meshes'
write(*,*) '        that contain triangles. THIS MODE IS VERY TIME-CONSUMING!'
write(*,*) ' '
write(*,*) '    se[t] [no]T0c[ompute]'
write(*,*) '        Controls whether initial temperatures calculation for ITMODE=3'
write(*,*) '        needs to be done or not. Facilitates output of parameters that'
write(*,*) '        where used to compute initial temperatures.'
write(*,*) ' '
write(*,*) '    se[t] [no]ve[rbose]'
write(*,*) '        Output mode. Typical in use when monitoring set.'
write(*,*) ' '
! write(*,*) '    se[t] ve[locity] value1 [unit1 = value2 unit2]'
! write(*,*) '        Functionality of this command depends on whether output is for GMT'
! write(*,*) '        or for P. For GMT, it defines a multiplication factor on velocity'
! write(*,*) '        data. For P, this command defines the scaling of velocities'
! write(*,*) '        to grid units. Unit fields are ignored on input,'
! write(*,*) '        but are convenient to keep track of the meaning of numbers. Example:'
! write(*,*) '                set velocity 1e-9 m/s = 5000 m'
! write(*,*) '        The first argument is also used to define the size of a scale bar.'
! write(*,*) '        (see also: `set gmt'', `set [no]scaling'') '
! write(*,*) ' '
write(*,*) '    se[t] visc[ale] factor'
write(*,*) '        Set scaling factor for calculating linearized viscosities'
write(*,*) '        (e.g. for transforming MPa.s to Pa.s)'
 else if (ofwhat(1:4) .eq. 'sedi') then
write(*,*) 'sedi[ment] [no[save]] file'
write(*,*) '    Write erosion surface coordinates to ascii output file. Argument'
write(*,*) '    [no]save is used to control if (saved) erosion surfaces or the'
write(*,*) '    current erosion surface needs to be written.'
write(*,*) '    (See also: `save'', `monitor'', `set verbose'') '
 else if (ofwhat(1:2) .eq. 'sh') then
write(*,*) 'sh[ow] parameter [> file]'
write(*,*) '    sh[ow] ti[me][(itime)]'
write(*,*) '        displays FE-output times'
write(*,*) '    sh[ow] calc'
write(*,*) '        displays times of thermal calculations'
write(*,*) '    sh[ow] fe'
write(*,*) '        shows finite-element characteristic numbers'
write(*,*) '    sh[ow] pl'
write(*,*) '        shows characteristic variables of program'
write(*,*) '    sh[ow] bc'
write(*,*) '        shows timing of various boundary conditions'
write(*,*) '    sh[ow] dt'
write(*,*) '        shows time step group characteristics'
 else if (ofwhat(1:5) .eq. 'slip ' .or. ofwhat(1:5).eq.'slipp') then
write(*,*) 'slip[pery] file'
write(*,*) '   Write slippery node directions of relative displacment to file.'
write(*,*) '   In gmt mode no symbol setting and scaling is needed.'
write(*,*) '   If set nogmt is on (plotting in P): Environmental variable PU must be set to'
write(*,*) '   the "P" subdirectory putil (e.g. setenv PU /local/P/putil).  (See also: `set'
write(*,*) '   gmt'', `set sysi'',`set deformed'').'
 else if (ofwhat(1:5) .eq. 'slips') then
write(*,*) 'slipsides filename'
write(*,*) '   Write slippery node weights as element contour data to a gmt file.'
write(*,*) '   See also: `conn''.' 
 else if (ofwhat(1:4) .eq. 'stra') then
write(*,*) 'stra[in] file [no[dev]]'
write(*,*) '    Write file for plotting element strains as arrows in principal stress'
write(*,*) '    directions using psvelo -Sx. (See also: `set [no]gmt'', `set ascii'', `set'
write(*,*) '    deformed'', `set elmi'', `monitor'', `set arrowhead'', `set verbose'', `set'
write(*,*) '    strain'') .'
 else if (ofwhat(1:4) .eq. 'stre') then
write(*,*) 'stre[ss] file [ [no]dev[iatoric] ] [tecton]'
write(*,*) '    Write file for plotting element stresse tensors as eigenvectors (x, y, sigma1, sigma2, angle)'
write(*,*) '    using GMT "psvelo -Sx" or "velo -Sx". Flags ""nodev"" and ""dev"" can be used to'
write(*,*) '    get non-deviatoric or deviatoric stresses. In case that the "tecton" flag'
write(*,*) '    is given, stresses are written as ASCII data in GTECTON pre-stress format.'
write(*,*) '    In "tecton" mode, giving the ""nodev"" argument will have the effect of'
write(*,*) '    outputing isotropic stress components only. (See also: "set [nogmt]", "set'
write(*,*) '    ascii", "set deformed", "set elmi", "monitor", "set arrowhead", "set'
write(*,*) '    verbose", "set stress" "set [no]scaling").'
 else if (ofwhat(1:4) .eq. 'surf') then
write(*,*) 'surf[lux] iside [file]'
write(*,*) '    Write surface heat flux magnitude data. First select elements'
write(*,*) 'with ''monitor elements''. ""iside"" is the side of the elements on'
write(*,*) 'which the surface heat flux is calculated.'
 else if (ofwhat(1:2) .eq. 'tc') then
write(*,*) 'tc file'
write(*,*) '    Write temperature in center of elements as contour data.'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'', `set cst'').'
 else if (ofwhat(1:2) .eq. 'te') then
write(*,*) 'te[mperatures] [tecton] file'
write(*,*) '    Writes temperature contour data to file. With argument ""tecton"" a tecton'
write(*,*) '    initial temperatures file is written. (See also: `set ascii'', `set'
write(*,*) '    deformed'', `sample'', `monitor'', `set npi'', `set verbose'', `set'
write(*,*) '    temperature'').'
 else if (ofwhat(1:1) .eq. 'u') then
write(*,*) 'u [number]'
write(*,*) '    undo save of last command(s).'
 else if (ofwhat(1:2) .eq. 've') then
write(*,*) 've[locity] [normal|slip|split|total] file'
write(*,*) '    Writes nodal velocity arrows to file.'
write(*,*) '    (See also: `sample'', `set ascii'', `set deformed'', `set npi'', `monitor'','
write(*,*) '    `set arrowhead'', `set verbose'', `set velocity'', `set [no]scaling'').'
 else if (ofwhat(1:2) .eq. 'vi') then
write(*,*) 'vi[scosity] file [no[log]] [no[nu[mbers]]] [no[li[near]] [no[raw]]'
write(*,*) '    writes viscosity contour data to file.'
write(*,*) '        log:  take 10-log of viscosity'
write(*,*) '      nolog:  no logaritmic scaling of viscosity'
write(*,*) '         nu:  write viscosity numbers in the form of ''pt''-commands to file'
write(*,*) '       nonu:  write viscosities as a function of spatial variables to file'
write(*,*) '         li:  calculate Newtonian equivalent of (non-linear) viscosity.'
write(*,*) '       noli:  no scaling to Newtonian viscosity'
write(*,*) '        raw:  original viscosities, without viscosity limit'
write(*,*) '      noraw:  processed viscosities, including viscosity limit'
write(*,*) '    There is no default; the last setting is used.'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'', `set cst'') .'
 else if (ofwhat(1:3) .eq. 'vol') then
write(*,*) 'vol[ume_change] [rate] file'
write(*,*) '    Write volume changes or volume change rates as contour data. (See also:'
write(*,*) '    `set ascii'', `set deformed'', `set triang'', `monitor'', `set verbose'','
write(*,*) '    `set cst'') .'
 else if (ofwhat(1:1) .eq. 'w' .and. len(trim(ofwhat)).eq.1) then
write(*,*) '?w says, lengths: ', len(ofwhat), len(trim(ofwhat))
write(*,*) 'w[!] [file]'
write(*,*) '    Write command buffer to command file. If ""w!"" do not overwrite original.'
 else if (ofwhat(1:4) .eq. 'wait') then
write(*,*) 'wait'
write(*,*) '    resumes execution until input from keyboard has been given.'
 else if (ofwhat(1:2) .eq. 'we') then
write(*,*) 'we [file]'
write(*,*) '    Write command buffer to command file and edit it'
 else if (ofwhat(1:2) .eq. 'wf') then
write(*,*) 'wf[emi] file'
write(*,*) '     Writes a f3d input file.'
write(*,*) '    (See also:) '
 else if (ofwhat(1:2) .eq. 'wq') then
write(*,*) 'wq [file]'
write(*,*) '    Write command buffer to command file and exit program'
 else if (ofwhat(1:2) .eq. 'xl') then
write(*,*) 'xl length'
write(*,*) '        Specifies horizontal plot dimension for geographic projection.'
write(*,*) '        (See also: ""set project"", ""gp"" ""lon"", ""lat"", ""yl"", ""pr"") .'
 else if (ofwhat(1:2) .eq. 'yi') then
write(*,*) 'yi[eld] file'
write(*,*) '    Write yield stress as contour data. (See also: `set ascii'','
write(*,*) '    `set deformed'', `set triang'', `monitor'', `set verbose'', `set cst'').'
 else if (ofwhat(1:2) .eq. 'yl') then
write(*,*) 'yl length'
write(*,*) '        Specifies vertical plot dimension for geographic projection.'
write(*,*) '        (See also: ""set project"", ""gp"" ""lon"", ""lat"", ""xl"", ""pr"").'
 else
write(*,*) 'argument not recognised. No manual entry for: ', trim(ofwhat), '.'
 endif

end subroutine

subroutine plt3d_manual(helpOfThis)

implicit none

integer :: nChars, iHelp

character(len=*) helpOfThis
character(len=20) ofwhat

integer :: helpListWidth, nHelpItems
#ifdef EXPERIMENTALTRANSIENTDISCREEP
character(len=14) :: helpList(48)
#else
character(len=14) :: helpList(47)
#endif

! The help request can be shorter than the
! the longer options
! ? help mo
! with an argument of length 2 should not crash on
! if (ofWhat(1:5).eq.estre) then

ofwhat = "                    "  ! twenty spaces

nChars           = len(helpOfThis)
ofwhat(1:nChars) = helpOfThis(1:nChars)

helpList = (/"             !",&
             "             #",&
             "             ?",&
             "            cf",&
             "          conn",&
             "           cou",&
             "            cp",&
             "           deb",&
             "            di",&
             "             e",&
             "           ebc",&
             "            ec",&
             "          elmc",&
             "           era",&
             "         error",&
             "         estra",&
             "         estre",&
             "          femi",&
             "          femo",&
             "           gri",&
             "             h",&
             "            ma",&
             "            mo",&
             "           nbc",&
             "          nomo",&
             "            np",&
             "             q",&
             "             r",&
             "            ra",&
             "           rot",&
             "            rt",&
             "          samp",&
             "            se",&
             "            sh",&
             "          stra",&
             "          stre",&
#ifdef EXPERIMENTALTRANSIENTDISCREEP
             "  taylorstress",&
#endif
             "            tc",&
             "            te",&
             "      traction",&
             "             u",&
             "            ve",&
             "            vi",&
             "             w",&
             "          wait",&
             "         write",&
             "            we",&
             "            wq"/)

if (ofwhat(1:4) .eq. 'list') then
    nHelpItems =  size(helpList)
    helpListWidth=4

    write(*,*) 'all available commands:'
    do iHelp=1,nHelpItems-1
        if (mod(iHelp,helpListWidth) .eq. 0) then
            write(*,'(a)'             ) helpList(iHelp)
        else
            write(*,'(a)',advance="no") helpList(iHelp)
        endif
    enddo

    write(*,'(a)'             ) helpList(iHelp)
    write(*,*) ' see ""? [command]"" for help on individual commands'
else if (ofwhat(1:1) .eq. '!') then
write(*,*) '! command'
write(*,*) '    Issue a shell command'
 else if (ofwhat(1:1) .eq. '#') then
write(*,*) '# comment'
write(*,*) '  records beginning with `#'' will be ignored.'
write(*,*) '  Inline comments start with `\"''. The rest of the record will be ignored.'
 else if (ofwhat(1:1) .eq. '?') then
write(*,*) '? [command]'
write(*,*) '        command,        gives verbose info on "command".'
 else if (ofwhat(1:2) .eq. 'cf') then
write(*,*) 'cf file'
write(*,*) '    Execute commands from command file `file''. (See also: `restart'').'
 else if (ofwhat(1:4) .eq. 'conn') then
write(*,*) 'conn[ections] file'
write(*,*) '    Write connections file in DX format (see also: "monitor", "set").'
 else if (ofwhat(1:3) .eq. 'cou') then
write(*,*) 'cou[lomb] [mhu] file'
write(*,*) '    writes vector data to file, specifying the normal vector of the optimally'
write(*,*) '    oriented fault plane. The magnitude of this vector is the associated shear'
write(*,*) '    stress. Coulomb stresses are calculated from the stress tensor and the'
write(*,*) '    friction coefficient "mhu".'
write(*,*) '    (See also: `set ascii'', `set deformed'', `monitor'', `set verbose'').'
 else if (ofwhat(1:2) .eq. 'cp') then
write(*,*) 'cp FromFile ToFile'
write(*,*) '    copies Fe-disk file FromFile to ToFile, converting binary to ascii'
write(*,*) '    or vice versa.'
 else if (ofwhat(1:3) .eq. 'deb') then
write(*,*) 'deb[orah] file'
write(*,*) '    writes Deborah number data to file.'
write(*,*) '    (See also: `set ascii'', `set deformed'', `monitor'', `set verbose'').'
 else if (ofwhat(1:2) .eq. 'di') then
write(*,*) 'di[splacements] file [total|normal|slippery|faulted]'
write(*,*) '     Writes nodal displacement vector data to file'
write(*,*) '     Selection of total, continuous, slippery or faulted node displacement is '
write(*,*) '     accomplished by giving the corresponding flag. (See also: "set ascii",'
write(*,*) '     "set deformed", "set np", "set npi", "monitor", "set verbose").'
 else if ((ofwhat(1:1) .eq. 'e' .and. len(trim(ofwhat)).eq.2) .or. &
           ofwhat(1:3) .eq.'edi') then
write(*,*) 'e[!] [file]'
write(*,*) '    edit command file. If "e!" use the original file.'
 else if (ofwhat(1:3) .eq. 'ebc') then
write(*,*) 'ebc w|s file'
write(*,*) '    Write element based mechanical boundary conditions to file. Flags:'
write(*,*) '    w: Winkler restoring pressure'
write(*,*) '    s: Stresses (outputs resolved tractions)'
write(*,*) '    (See also: `nbc'', `set sysi'', `set deformed'', `set noscaling'').'
 else if (ofwhat(1:2) .eq. 'ec') then
write(*,*) 'ec[ho] [-n] [text]'
write(*,*) '    echoes text. Flag -n suppresses line feed.'
 else if (ofwhat(1:4) .eq. 'elmc') then
write(*,*) 'elmc[enter] file'
write(*,*) '    Write element center coordinates to file. (See also: `set ascii'','
write(*,*) '    `set deformed'', `set elmi'', `monitor'', `set verbose'').'
 else if (ofwhat(1:3) .eq. 'era') then
write(*,*) 'era[te] [flag] file'
write(*,*) '    Write effective total strain rates or total strain rate components as'
write(*,*) '    contour data. Produces a DX file which includes the nodal point'
write(*,*) '    coordinate file and the connectivity file. Therefore, the "np" and "conn"'
write(*,*) '    commands should be issued to generate these include files, before writing'
write(*,*) '    effective rates to file.'
write(*,*) ' '
write(*,*) '    Possible flags include:'
write(*,*) '        e       effective values'
write(*,*) '        xx      xx-component of deviator'
write(*,*) '        yy      yy-component of deviator'
write(*,*) '        zz      zz-component of deviator'
write(*,*) '        xy      xy component'
write(*,*) '        xz      xz component'
write(*,*) '        yz      yz component'
write(*,*) ' '
write(*,*) '    (See also: "set ascii", "set deformed", "np", "conn", "monitor",'
write(*,*) '    "set verbose").'
 else if (ofwhat(1:5) .eq. 'error') then
write(*,*) 'error file'
write(*,*) '    Writes the norm of the a posteriori error vector'
 else if (ofwhat(1:5) .eq. 'estra') then
write(*,*) 'estra[in] [flag] file'
write(*,*) '    Write effective strain or deviatoric strain components as contour data.'
write(*,*) '    Produces a DX file which includes the nodal point coordinate file and the'
write(*,*) '    connectivity file. Therefore, the "np" and "conn" commands should be'
write(*,*) '    issued to generate these include files, before writing effective strains'
write(*,*) '    to file.'
write(*,*) ' '
write(*,*) '    Possible flags include:'
write(*,*) '        e       effective values'
write(*,*) '        xx      xx-component of deviator'
write(*,*) '        yy      yy-component of deviator'
write(*,*) '        zz      zz-component of deviator'
write(*,*) '        xy      xy component'
write(*,*) '        xz      xz component'
write(*,*) '        yz      yz component'
write(*,*) ' '
write(*,*) '    (See also: "set ascii", "set deformed", "np", "conn", "monitor",'
write(*,*) '    "set verbose").'
 else if (ofwhat(1:5) .eq. 'estre') then
write(*,*) 'estre[ss] [flag] file'
write(*,*) '    Write effective stress or deviatoric stress components as contour data.'
write(*,*) '    Produces a DX file which includes the nodal point coordinate file and the'
write(*,*) '    connectivity file. Therefore, the "np" and "conn" commands should be'
write(*,*) '    issued to generate these include files, before writing effective stresses'
write(*,*) '    to file.'
write(*,*) '  '
write(*,*) '    Possible flags include:'
write(*,*) '        e       effective values'
write(*,*) '        xx      xx-component of deviator'
write(*,*) '        yy      yy-component of deviator'
write(*,*) '        zz      zz-component of deviator'
write(*,*) '        xy      xy component'
write(*,*) '        xz      xz component'
write(*,*) '        yz      yz component'
write(*,*) '  '
write(*,*) '    (See also: "set ascii", "set deformed", "np", "conn", "monitor", "set'
write(*,*) '    verbose").'
 else if (ofwhat(1:4) .eq. 'femi') then
write(*,*) 'femi[n] [MechFile [TempFile]]'
write(*,*) '    Read finite element (TECTON) input file(s). "MechFile" (typically named'
write(*,*) '    "TECIN.DAT") contains mesh geometry, material properties and boundary'
write(*,*) '    conditions for the mechanical problem. If neccesary, "TempFile" (typically'
write(*,*) '    named "TEMPIN.DAT") is read. "TempFile" contains material properties and'
write(*,*) '    boundary conditions relevant for the thermal problem. (See also: `femo'').'
 else if (ofwhat(1:4) .eq. 'femo') then
write(*,*) 'femo[ut] [file] [bin[ary]|asc[ii]]'
write(*,*) '    Read finite element output files. "file" is the generic finite element output'
write(*,*) '    file which should be similar to the output file specified for pln being of the'
write(*,*) '    form "[dirname/][basename]". If nothing is specified the default "FEDSK.DAT" is'
write(*,*) '    used. Normally the input format (binary or ascii) is determined automatically.'
write(*,*) '    The file format can be forced by giving "binary" or "ascii" arguments.'
write(*,*) '    (See also: `femi'').'
 else if (ofwhat(1:2) .eq. 'fl') then
write(*,*) 'fl[ux] file'
write(*,*) '    Write element heat fluxes to a DX file. (See also: `set ascii'','
write(*,*) '    `set deformed'', `set elmi'', `monitor'', `set arrowhead'', `set verbose'','
write(*,*) '    `set flux'').'
 else if (ofwhat(1:2) .eq. 'gr') then
write(*,*) 'gr[id] file'
write(*,*) '    Generates a "DX" datafile of the finite element mesh. The method is to'
write(*,*) '    include nodal point data, written to disk using the "np" command,'
write(*,*) '    and connectivity data, written using the "conn" command into a'
write(*,*) '    DataExplorer file. It is, therefore, neccessary to first perform the'
write(*,*) '    "np file" and "conn file" commands, before assembly of grid data can'
write(*,*) '    take place (See also, "np", "conn").'
 else if (ofwhat(1:1) .eq. 'h') then
write(*,*) 'h[istory] [number]'
write(*,*) '    prints history of last commands.'
 else if (ofwhat(1:2) .eq. 'ma') then
write(*,*) 'ma[terial] file'
write(*,*) '    Write material numbers as contour data. Produces a DX file which includes'
write(*,*) '    the nodal point coordinate file and the connectivity file. Therefore, the'
write(*,*) '    "np" and "conn" commands should be issued to generate these include files,'
write(*,*) '    before writing material numbers to file.'
write(*,*) ' '
write(*,*) '    (See also: "set ascii", "set deformed", "np", "conn", "monitor", "set'
write(*,*) '    verbose").'
 else if (ofwhat(1:2) .eq. 'mo') then
write(*,*) 'mo[nitor] flag number[s]'
write(*,*) '    monitor output at a specific node or element number. There is no'
write(*,*) '    default. Turned off by `nomonitor nodes'' or `nomonitor elements''.'
write(*,*) '        This command affects the syntax of commands that normally require a'
write(*,*) '    file name argument. If, in these commands, the file name argument is'
write(*,*) '    omitted, information is written to standard output. If the file name'
write(*,*) '    argument is given, monitor data are written (ASCII format) to this'
write(*,*) '    file. If `noverbose'' has been set, monitor information is very compact.'
write(*,*) '        "flag" should be one of:'
write(*,*) '            n[odes],            NUMBERS are interpreted as nodal point no''s.'
write(*,*) '            e[lements],         NUMBERS are interpreted as element no''s.'
write(*,*) '    Ranges of numbers can be given as "istart,istep,iend" or "istart,iend"'
write(*,*) '    A valid statement is therefore:'
write(*,*) '        monitor e 3,5,28'
write(*,*) '  '
write(*,*) '    (See also: `set verbose'').'
 else if (ofwhat(1:3) .eq. 'nbc') then
write(*,*) 'nbc [a|d|v|f|i|r] file'
write(*,*) '    Write mechanical nodal boundary condition to file. Flags:'
write(*,*) '    a: all boundary condition types'
write(*,*) '    d: displacement boundary condition'
write(*,*) '    v: velocity boundary condition'
write(*,*) '    f: force boundary condition'
write(*,*) '    i: initial displacement boundary condition'
write(*,*) '    r: strain rate boundary condition'
write(*,*) '    (See also: `set ascii'', `set append'', `set sysi'', `set deformed'','
write(*,*) '    `set noscaling'').'
 else if (ofwhat(1:4) .eq. 'nomo') then
write(*,*) 'nomo[nitor] n[ode]|e[lement]|r[ow]|c[olumn]'
write(*,*) '    Turn monitoring of nodes, elements, rows or columns off.'
 else if (ofwhat(1:2) .eq. 'np') then
write(*,*) 'np file'
write(*,*) '    Write nodal coordinates to a DataExplorer file. (See also: "conn",'
write(*,*) '    "grid", "monitor", "set deformed", "set verbose").'
 else if (ofwhat(1:1) .eq. 'q') then
write(*,*) 'q[uit]'
write(*,*) '    stops execution of program'
 else if (ofwhat(1:1) .eq. 'r' .and. len(trim(ofwhat)).eq.1) then
write(*,*) 'r'
write(*,*) '    restart program. If "r!" restart from original command file. Most variable'
write(*,*) '    settings are remembered. (See also: `cf'').'
 else if (ofwhat(1:2) .eq. 'ra') then
write(*,*) 'ra[te] file'
write(*,*) '    Write strain rate tensors to a DX data file. Produces a DX file which'
write(*,*) '    includes the nodal point coordinate file and the connectivity file.'
write(*,*) '    Therefore, the "np" and "conn" commands should be issued to generate these'
write(*,*) '    include files, before writing strain rate tensors to file.  (See also: "set'
write(*,*) '    ascii", "set deformed", "np", "conn", "monitor", "set verbose", "set'
write(*,*) '    strain").'
 else if (ofwhat(1:3) .eq. 'rot') then
write(*,*) 'rot[ation] [no[rate] file'
write(*,*) '     Writes (infinitesimal) rotation tensors or rotations rate tensor to file.'
write(*,*) '     (See also: `set ascii'', `set deformed'', `monitor'', `set verbose'').'
 else if (ofwhat(1:2) .eq. 'rt') then
write(*,*) 'rt [stress|strain|rate] [side#]'
write(*,*) '    Echoes normal and shear components of tensors as resolved on faces of'
write(*,*) '    elements. Only works with monitor elements on ("mon e")'
write(*,*) '    (See also: `set deformed'', `monitor'', `set verbose'').'
 else if (ofwhat(1:4) .eq. 'samp') then
write(*,*) 'samp[le] d|v|T|s|e|r|c|f coordinate [filename]'
write(*,*) '   compute quantities at coordinate[s]'
write(*,*) '   d = displacement'
write(*,*) '   v = velocity'
write(*,*) '   T = temperature'
write(*,*) '   s = stress'
write(*,*) '   e = strain'
write(*,*) '   r = strain rate'
write(*,*) '   c = viscosity'
write(*,*) '   f = heat flux'
write(*,*) '   A single coordinate is specified as x y z. A range of coordinates is specified'
write(*,*) '   as xmin,xmax[,xinc] ymin,ymax[,yinc] zmin,zmax[,zinc].'
write(*,*) '   Sample is specifically used to interpolate (based on shape functions) nodal'
write(*,*) '   quantities. Use "monitor np" in combination with "disp", "vel" or "temp" to'
write(*,*) '   get nodal values. For element quantities, functionality of sample partly'
write(*,*) '   overlaps with that of "monitor e" in combination with "stress", "strain",'
write(*,*) '   "rate" or "flux".'
 else if ((ofwhat(1:3) .eq. 'se ').or.(ofwhat .eq. 'set')) then
write(*,*) 'se[t] parameter [value[s]]'
write(*,*) '  '
write(*,*) '    se[t] ar[rowhead] size angle [no[co[nt]]]'
write(*,*) '        Define size and top angle of arrows.'
write(*,*) ' '
write(*,*) '    se[t] [no]ap[pend]'
write(*,*) '        Append data to file.'
write(*,*) ' '
write(*,*) '    se[t] [no]as[cii]'
write(*,*) '        Output file format'
write(*,*) ' '
! write(*,*) '    se[t] bcsy[mbol] type number'
! write(*,*) '        Set symbol number of boundary condition of type `type''.'
! write(*,*) '        `type'' can be one of;'
! write(*,*) '        fi[xed]:        for displacements fixed in more than one direction'
! write(*,*) '        sl[ide]:        for displacements fixed in one direction'
! write(*,*) '        fo[rce]:        for nodal force boundary conditions'
! write(*,*) '        ve[locity]:     for nodal velocity boundary conditions'
! write(*,*) '        st[ress]:       for element stress boundary conditions'
! write(*,*) '        wi[nkler]:      for nodal Winkler boundary conditions'
! write(*,*) ' '
write(*,*) '    se[t] [no]bi[nary]'
write(*,*) '        Output file format'
write(*,*) ' '
write(*,*) '    se[t] bou[guer] [Yref]'
write(*,*) '        Sets gravity calculation mode to "bouguer". In this mode, the'
write(*,*) '        gravity attraction of masses below the horizon Y=Yref are calculated.'
write(*,*) '        Yref is set to 0 initially, and may be reset by both "set bouguer"'
write(*,*) '        and "set free-air" commands. (See also "set raw", set "free-air",'
write(*,*) '        "gref", "grav", "set expand", "set density").'
write(*,*) ' '
write(*,*) '    se[t] [no]de[formed] [amplification]'
write(*,*) '       Original or deformed mesh output. Amplification factor can be used'
write(*,*) '        to exaggerate deformation.'
write(*,*) '  '
write(*,*) '    se[t] den[sity] number'
write(*,*) '        Sets density of elements which are y-expanded to a depth of 300 km.'
write(*,*) '        This option is used for calculating gravity anomalies. If the density'
write(*,*) '        is set to a value less or equal zero, current elements densities'
write(*,*) '        will be used (See also "set raw", set "free-air", "set bouguer",'
write(*,*) '        "gref", "grav", "set expand").'
write(*,*) ' '
! write(*,*) '    se[t] di[splacement] value1 unit1 = value2 unit2'
! write(*,*) '        Define displacements scaling to grid units.'
! write(*,*) ' '
write(*,*) '    se[t] [no]ec[ho] [number]'
write(*,*) '        Set echoing mode. `noecho'' is identical to `echo 0''.'
write(*,*) '        echo 1:   verbose mode.'
write(*,*) '        echo 2:   input files reading.'
write(*,*) '        echo 3,4: program calculations'
write(*,*) '        echo 5:   verbose program calculations'
write(*,*) ' '
write(*,*) '    se[t] elmi[nterval] [base] number'
write(*,*) '        Sets element base and increment for plotting of arrows.'
write(*,*) '  '
write(*,*) '    se[t] ex[pand] x|y iside numbers'
write(*,*) '        Sets element numbers and element side numbers which have to be'
write(*,*) '        expanded horizontally or vertically during gravity calculations.'
write(*,*) '        Expansion is done by 10 times maximum x- or y-values.'
write(*,*) '        (See also "gref", "grav", "set raw", "set bouguer", "set free-air").'
write(*,*) ' '
! write(*,*) '    se[t] fl[ux] value1 unit1 = value2 unit2'
! write(*,*) '        Define heat flux scaling to grid units.'
! write(*,*) ' '
write(*,*) '    se[t] fre[e-air] [Yref]'
write(*,*) '        Sets gravity calculation mode to "free-air". In this mode,'
write(*,*) '        the gravity attraction of masses below the surface are calculated.'
write(*,*) '        Subsequently, the signal is reduced to the horizon Y=Yref.'
write(*,*) '        Yref is set to 0 initially, and may be reset by both "set bouguer"'
write(*,*) '        and "set free-air" commands. (See also "set bouguer", "set raw",'
write(*,*) '        "gref", "grav", "set expand", "set density").'
write(*,*) ' '
write(*,*) '    se[t] it[ime] number '
write(*,*) '        Set output time step number. Increment "++" and decrement "--"'
write(*,*) '        operators are also allowed.'
write(*,*) '  '
! write(*,*) '    se[t] np n|e'
! write(*,*) '        Sets flag to compute/output nodal quantities on the basis of nodes'
! write(*,*) '        (n) or elements (e). This can be important in the case where slippery'
! write(*,*) '        nodes or faulted nodes occur.'
! write(*,*) ' '
write(*,*) '    se[t] npi[nterval] [base] number'
write(*,*) '        Sets nodal point base and increment for plotting of arrows.'
write(*,*) ' '
write(*,*) '    se[t] opn [flag [NUMZ [scale [axnu]]]]'
write(*,*) '        Sets "out-of-plane" mode. These settings define the way opn fe-input'
write(*,*) '        files are converted to f3d files. "flag" determines whether grid'
write(*,*) '        dimensions ("x") or displacements ("d") are used to specify out of'
write(*,*) '        plane dimensions. NUMNZ specifies the number of additional nodal point'
write(*,*) '        planes which is being generated. "scale" is a linear scaling parameter.'
write(*,*) '        "axnu" specifies the fraction of the derived dimensions that is being'
write(*,*) '        added to the maximum and subtracted from the minimum dimension values.'
write(*,*) '        Example: set opn d 4 1.0 0.05'
write(*,*) ' '
! write(*,*) '    se[t] ra[te] value1 unit1 = value2 unit2'
! write(*,*) '        Define strain rate scaling to grid units.'
! write(*,*) ' '
! write(*,*) '    se[t] [no]sc[aling]'
! write(*,*) '        Turns boundary condition symbols automatic scaling on or off.'
! write(*,*) ' '
write(*,*) '    se[t] raw'
write(*,*) '        Sets gravity calculation mode to "raw". Gravity calculations, performed'
write(*,*) '        in this mode will result in the gravity anomaly as if measured along'
write(*,*) '        the current model surface. (See also "set bouguer", "set free-air",'
write(*,*) '        "grav", "gref", "set expand", "set density").'
write(*,*) ' '
write(*,*) '    se[t] sh[rink] [factor]'
write(*,*) '        Shrinks elements by a factor. Useful for testing nodal point'
write(*,*) '        connectivity. (See also "grid").'
write(*,*) ' '
! write(*,*) '    se[t] stra[in] value unit1 = value2 unit2'
! write(*,*) '        Define strain scaling to grid units.'
! write(*,*) ' '
! write(*,*) '    se[t] stre[ss] value comment = unit comment'
! write(*,*) '        Define stress scaling to grid units.'
! write(*,*) '        "comment"-fields are ignored on input, but are convenient to keep'
! write(*,*) '        track of the meaning of numbers. Example:'
! write(*,*) '                set stress 2.E8 Pa = 5.E3 m'
! write(*,*) '        The first argument is also used to define the size of a scale bar.'
! write(*,*) ' '
! write(*,*) '    se[t] sy[mbolsize] size'
! write(*,*) '        Set symbol size in plot units.'
! write(*,*) ' '
write(*,*) '    se[t] tag [number]'
write(*,*) '        Set the material tag of the subdomain which is'
write(*,*) '        to be written to file. Formats are:'
write(*,*) '        displacement: idx, x, y, z, dx, dy, dz'
write(*,*) '        velocity:     idx, x, y, z, vx, vy, vz'
write(*,*) '        traction:     idx, nx, ny, ny, tx, ty, tz'
write(*,*) '                    with nx, ny, nz the normal vector to the triangl'
write(*,*) ' '
write(*,*) '    se[t] tem[perature] a[bsolute]|r[elative]|i[ncrement] [itime]'
write(*,*) '        Set temperatures mode.'
write(*,*) '            `absolute''  =current temperatures'
write(*,*) '            `relative''  =temperature differences relative to some other'
write(*,*) '                         time, specified by itime (default: itime=0)'
write(*,*) '            `increment'' =temperature differences relative to previous output'
write(*,*) ' '
write(*,*) '    se[t] [no]test'
write(*,*) '        Perform additional tests on GTECTON input data'
write(*,*) '        Default: notest'
write(*,*) ' '
write(*,*) '    se[t] [no]tr[iangulate]'
write(*,*) '        Triangulation mode. Usefull when contouring on meshes'
write(*,*) '        that contain triangles. THIS MODE IS VERY TIME-CONSUMING!'
write(*,*) ' '
write(*,*) '    se[t] [no]T0c[ompute]'
write(*,*) '        Controls whether initial temperatures calculation for ITMODE=3'
write(*,*) '        needs to be done or not. Facilitates output of parameters that'
write(*,*) '        where used to compute initial temperatures.'
write(*,*) ' '
write(*,*) '    se[t] [no]ve[rbose]'
write(*,*) '        Output mode. Typical in use when monitoring set.'
write(*,*) ' '
write(*,*) '    se[t] nosoelim or dosoelim'
write(*,*) ' '
write(*,*) '    se[t] visc[ale] factor'
write(*,*) '        Set scaling factor for calculating linearized viscosities'
write(*,*) '        (e.g. for transforming MPa.s to Pa.s)'
write(*,*) ' '
write(*,*) '    se[t] vtu [option] [on/off]'
write(*,*) '        Option only relevant for writing output of 3D runs.'
write(*,*) '        To prevent the output from becoming overly large,'
write(*,*) '        the blocks of content can be turned on or off.'
write(*,*) '        Possible options are:'
write(*,*) ' '
write(*,*) '        - prestress;        default: ''off'' '
write(*,*) '        - materials;        default: ''off'' '
write(*,*) '        - velocity;         default: ''on''  '
write(*,*) '        - displacement;     default: ''on''  '
write(*,*) '        - stress;           default: ''off'' '
write(*,*) '        - strain;           default: ''off'' '
#ifdef EXPERIMENTALTRANSIENTDISCREEP
write(*,*) '        - taylorstress;     default: ''off'' '
#endif
write(*,*) '        - strainrate;       default: ''off'' '
write(*,*) '        - elempartition;    default: ''off'' '
write(*,*) '        - vertexpartition;  default: ''off'' '
write(*,*) '        - writeASCII;       [still in test phase. fixed in ''on'']'
write(*,*) '        - elemQualitySide   default: ''off'' '
write(*,*) '        - elemQualityFace   default: ''off'' '
write(*,*) '        - elemQualityVolume default: ''off'' '
write(*,*) '        - velocityBC        default: ''off'' '
write(*,*) '        - slideInPlainBC    default: ''off'' '
write(*,*) '        - slideOnLineBC     default: ''off'' '
write(*,*) '        - fixedPointBC      default: ''off'' '
write(*,*) '        - displacementBC    default: ''off'' '
write(*,*) '        - forceBC           default: ''off'' '
write(*,*) '        - EulerAnglesBC     default: ''off'' '
write(*,*) '        - error             default: ''off'' '
write(*,*) ' '
write(*,*) '        Note that ''1'' and ''0'' are also accepted for ''on'' and ''off''.'
write(*,*) ' '
write(*,*) '        elemQualitySide   gives the length of the shortest side divided by'
write(*,*) '                          the length of the longest side.'
write(*,*) '        elemQualityFace   gives the area of the the smallest face divided by'
write(*,*) '                          the area of the largest face.'
write(*,*) '        elemQualityVolume gives the volume of the incribed sphere divided by'
write(*,*) '                          the volume of the circumscribed sphere.'
write(*,*) '        slideInPlainBC    gives arrows of unit length that are normal to the'
write(*,*) '                          plane in which the node is allowed to slide,' 
write(*,*) '                          in the positive direction.'
write(*,*) '        slideOnLineBC     gives arrows of unit length that are along the'
write(*,*) '                          line on which the node is allowed to slide,'
write(*,*) '                          in the positive direction.'
write(*,*) '        fixedPointBC      gives (1,1,1) arrows'
write(*,*) '        EulerAnglesBC     gives a system of global arrows '
write(*,*) '                          (3,0,0) (0,2,0) and (0,0,1)'
write(*,*) '                          rotated into local orientation.'
write(*,*) ' '
write(*,*) '        All options can be turned on or off simultaneously with the command:'
write(*,*) '        set vtu all on'
write(*,*) '        and'
write(*,*) '        set vtu all off'
write(*,*) '        respectively.'
 else if (ofwhat(1:2) .eq. 'sh') then
write(*,*) 'sh[ow] parameter [> file]'
write(*,*) '    sh[ow] ti[me][(itime)]'
write(*,*) '        displays FE-output times'
write(*,*) '    sh[ow] calc'
write(*,*) '        displays times of thermal calculations'
write(*,*) '    sh[ow] fe'
write(*,*) '        shows finite-element characteristic numbers'
write(*,*) '    sh[ow] pl'
write(*,*) '        shows characteristic variables of program'
write(*,*) '    sh[ow] bc'
write(*,*) '        shows timing of various boundary conditions'
write(*,*) '    sh[ow] dt'
write(*,*) '       shows time step group characteristics'
 else if (ofwhat(1:4) .eq. 'stra') then
write(*,*) 'stra[in] file'
write(*,*) '    Write strain tensors to a DX data file. Produces a DX file which includes'
write(*,*) '    the nodal point coordinate file and the connectivity file. Therefore,'
write(*,*) '    the "np" and "conn" commands should be issued to generate these include'
write(*,*) '    files, before writing strain tensors to file.  (See also: "set ascii", "set'
write(*,*) '    deformed", "np", "conn", "monitor", "set verbose", "set strain").'
 else if (ofwhat(1:4) .eq. 'stre') then
write(*,*) 'stre[ss] file [tecton]'
write(*,*) '    Write stress tensors to a DX data file. Produces a DX file which includes'
write(*,*) '    the nodal point coordinate file and the connectivity file. Therefore,'
write(*,*) '    the "np" and "conn" commands should be issued to generate these include'
write(*,*) '    files, before writing stress tensors to file. If "tecton" argument is given,'
write(*,*) '    stresses are written in TECTON pre-stress data format (See also: "set'
write(*,*) '    ascii", "set deformed", "np", "conn", "monitor", "set verbose", "set'
write(*,*) '    stress").'
#ifdef EXPERIMENTALTRANSIENTDISCREEP
 else if (ofwhat(1:12) .eq. 'taylorstress') then
write(*,*) 'taylorstress file [tecton]'
write(*,*) '    Write taylorstress [tensors] to a DX data file. Produces a DX file which includes'
write(*,*) '    the nodal point coordinate file and the connectivity file. Therefore,'
write(*,*) '    the "np" and "conn" commands should be issued to generate these include'
write(*,*) '    files, before writing taylorstress tensors to file. If "tecton" argument is given,'
write(*,*) '    taylorstresses are written in TECTON pre-stress data format (See also: "set'
write(*,*) '    ascii", "set deformed", "np", "conn", "monitor", "set verbose", "set'
write(*,*) '    taylorstress").'
#endif
 else if (ofwhat(1:2) .eq. 'tc') then
write(*,*) 'tc file'
write(*,*) '    Write temperature in center of elements as contour data.'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'','
write(*,*) '    `set verbose'').'
 else if (ofwhat(1:2) .eq. 'te') then
write(*,*) 'te[mperatures] [tecton] file'
write(*,*) '    Writes temperature contour data to file. Writing nodal points coordinates'
write(*,*) '    and connectivity is required before issuing this command. With argument'
write(*,*) '    "tecton" a tecton initial temperatures file is written. In this case,'
write(*,*) '    writing nodal coordinates and connectivity is not neccessary.'
write(*,*) '    (See also: `set ascii'', `set deformed'', `set triang'', `monitor'', `set npi'','
write(*,*) '    `set verbose'', `set temperature'').'
 else if (ofWhat(1:8) .eq. 'traction') then
write(*,*) 'traction file'
write(*,*) '    Writes traction and related quantities (stress, surface normals) to file. These'
write(*,*) '    quantities are meant for further processing before plotting.'
write(*,*) '    The quantities belong to element sides located on surfaces previously selected via'
write(*,*) '    set tag. Where two sets of elements share the same surfaces (i.e., where slippery'
write(*,*) '    or faulted nodes are present), quantities belonging to both are written.'
write(*,*) '    Each output line contains the following quantities:'
write(*,*) '    sxx syy szz sxy sxz syz x1 y1 z1 x2 y2 z2 x3 y3 z3 nx ny nz tx ty tz,'
write(*,*) '    where sij are the components of the total stress tensor, xi/yi/zi are the'
write(*,*) '    coordinates of the ith node of the element side, and ni and ti are the components'
write(*,*) '    ofthe surface-normal unit vector and of the traction vector, respectively.'
 else if (ofwhat(1:1) .eq. 'u') then
write(*,*) 'u [number]'
write(*,*) '    undo save of last command(s).'
 else if (ofwhat(1:2) .eq. 've') then
write(*,*) 've[locity] file [total|normal|slippery|faulted]'
write(*,*) '     Writes nodal velocities to file. Selection of total, continuous,'
write(*,*) '     slippery or faulted node velocity is accomplished by giving the'
write(*,*) '     corresponding flag. (See also: `set ascii'', `set deformed'', `set np'','
write(*,*) '     `set npi'', `monitor'', `set verbose'').'
 else if (ofwhat(1:2) .eq. 'vi') then
write(*,*) 'vi[scosity] [no[log]] [no[li[near]] [no[raw]] file'
write(*,*) '    writes viscosity contour data to file.'
write(*,*) '        log:  take 10-log of viscosity'
write(*,*) '      nolog:  no logaritmic scaling of viscosity'
write(*,*) '         li:  calculate Newtonian equivalent of (non-linear) viscosity.'
write(*,*) '       noli:  no scaling to Newtonian viscosity'
write(*,*) '        raw:  original viscosities, without viscosity limit'
write(*,*) '      noraw:  processed viscosities, including viscosity limit'
write(*,*) '    There is no default; the last setting is used.'
write(*,*) '    (See also: `set ascii'', `set deformed'', `monitor'', `set verbose'').'
 else if (ofwhat(1:8) .eq. 'vtuwrite') then
write(*,*) 'vtuwrite filename'
write(*,*) '    Writes the selected content from a 3D run to file in VTU format,'
write(*,*) '    which can be read and visualised by ParaView'
write(*,*) '    Content is selected using ''set vtu [option] [on/off]''.'
write(*,*) '    see: ''? set''.'
 else if (ofwhat(1:1) .eq. 'w' .and. len(trim(ofwhat)).eq.1) then
write(*,*) 'w[!] [file]'
write(*,*) '    Write command buffer to command file. If "w!" do not overwrite original.'
 else if (ofwhat(1:4) .eq. 'wait') then
write(*,*) 'wait'
write(*,*) '    resumes execution until input from keyboard has been given.'
 else if (ofwhat(1:2) .eq. 'we') then
write(*,*) 'we [file]'
write(*,*) '    Write command buffer to command file and edit it'
 else if (ofwhat(1:5) .eq. 'write') then
write(*,*) 'write'
write(*,*) '    Write file data.vtu that can be read by Paraview.'
write(*,*) '    See: `? set'' for `se[t] vtu'' on how to turn switches on and off'
write(*,*) '    to specify the content of this file.'
 else if (ofwhat(1:2) .eq. 'wq') then
write(*,*) 'wq [file]'
write(*,*) '    Write command buffer to command file and exit program'
 else
write(*,*) 'argument not recognised. No manual entry for: ', trim(ofwhat), '.'
 endif



end subroutine



end module
