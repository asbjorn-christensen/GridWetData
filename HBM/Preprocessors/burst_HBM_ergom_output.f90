        program burst_HBM_ergom_output
!========================================================================================================================
!       Burst HBM/ERGOM v5 output into constituents
!
!       Usage: burstHBM [p|b] <physfile|biofile> [prop*]
!
!       where 
!          p : means the file is a physics file (no name pattern parsing)
!          b : means the file is a biology file (no name pattern parsing)
!       and
!          prop ==  t|s|u|v|z for physics layers or z1|z2|ch  for biology layers (more to follow)
!
!       (for windows users the executable will be burstHBM.exe)  
!
!              ------------- physics layers -------------------------
!
!              t : extract temperature (3D)
!              s : extract salinity (3D)
!              u : extract meridional currents (3D)
!              v : extract zonal currents (3D)
!              z : extract sea surface elevation  (2D)
!
!              ------------- biology layers -------------------------
!
!              z1 : extract micro zooplankton (3D)
!              z2 : extract meso  zooplankton (3D)
!              ch : extract chlorophyll (3D)
!        
!       The current version does not consistently return non-zero for anormal execusions (but info printet to stdout is updated)  
!       ----------------------------------------------------------------------------------------------------------
!       Based on HBM reading example cmod2cdf.f90/bio2cdf.f90, provided by Zhenwen Wan [zw@dmi.dk] May 26, 2015
!       Reading hardcoded to current 4-area setup, writing customizable 
!
!       index  tag  domain
!         0     ns  North Sea 
!         1     id  Inner Danish Waters
!         2     ws  Wadden Sea 
!         3     bs  Baltic Sea 
!
!       The operational ERGOM v5 is dfferent from OPEC version. There are 12 three-dimensional state variables and 
!       2 sediment variables in the operational ERGOM v5.
!       
!       Now,  12 three-dimensional state variables refer to as
!                   nh4(i,j,k)=eco(1,id3d(jm+1-j,i,k))/scale               
!                   no3(i,j,k)=eco(2,id3d(jm+1-j,i,k))/scale
!                   po4(i,j,k)=eco(3,id3d(jm+1-j,i,k))/scale               
!                   Sio4(i,j,k)=eco(4,id3d(jm+1-j,i,k))/scale               
!                   dia(i,j,k)=eco(5,id3d(jm+1-j,i,k))/scale
!                   fla(i,j,k)=eco(6,id3d(jm+1-j,i,k))/scale
!                   cya(i,j,k)=eco(7,id3d(jm+1-j,i,k))/scale
!                   zo1(i,j,k)=eco(8,id3d(jm+1-j,i,k))/scale
!                   zo2(i,j,k)=eco(9,id3d(jm+1-j,i,k))/scale
!                   SPM1(i,j,k)=eco(10,id3d(jm+1-j,i,k))/scale ! suspended particulate matters in content of nitrogen
!                   SPM2(i,j,k)=eco(11,id3d(jm+1-j,i,k))/scale ! suspended particulate matters in content of silicate
!                   oxy(i,j,k)=eco(12,id3d(jm+1-j,i,k))/scale
!       
!       2 sediment variables refer to as
!                   SPM_sed1(i,j)=sed(1,id3d(jm+1-j,i,1))/scale               
!                   SPM_sed2(i,j)=sed(2,id3d(jm+1-j,i,1))/scale
!========================================================================================================================  
!       Linux build : 
!            ifort -o burstHBM -convert big_endian -e90 -I/usr/local/include -L/usr/local/lib -lnetcdff -lnetcdf  burst_HBM_ergom_output.f90
!
!       cross-compile to static windows executable: e.g.
!            NETCDF=/home/asbjorn/DTU/Ballastvand_SRA/IBMlib_port_to_windows/mingw_netcdf/working/NETCDF
!            /usr/bin/x86_64-w64-mingw32-gfortran -c -fconvert=big-endian -I${NETCDF}/include burst_HBM_ergom_output.f90
!            /usr/bin/x86_64-w64-mingw32-gfortran -static -o burstHBM.exe burst_HBM_ergom_output.o -L${NETCDF}/lib  -lnetcdf -lnetcdff -lnetcdf -lnetcdff
!========================================================================================================================
        use netcdf   
        implicit None
        
        integer,parameter :: n3d3nm0=479081,  n2d3nm0=18908
        integer,parameter :: n3d3nm1=1584802, n2d3nm1=80885
        integer,parameter :: n3d3nm2=103441,  n2d3nm2=11581
        integer,parameter :: n3d3nm3=6112717, n2d3nm3=119206
!       ---- physics ----
        integer(4),dimension(0:n3d3nm0),target :: uu0,vv0,ww0,ss0,tt0
        integer(4),dimension(0:n3d3nm1),target :: uu1,vv1,ww1,ss1,tt1
        integer(4),dimension(0:n3d3nm2),target :: uu2,vv2,ww2,ss2,tt2
        integer(4),dimension(0:n3d3nm3),target :: uu3,vv3,ww3,ss3,tt3
        integer(4),dimension(0:n2d3nm0),target :: wu0,wv0,zz0
        integer(4),dimension(0:n2d3nm1),target :: wu1,wv1,zz1
        integer(4),dimension(0:n2d3nm2),target :: wu2,wv2,zz2
        integer(4),dimension(0:n2d3nm3),target :: wu3,wv3,zz3
!       ---- ecovars ----
        integer,parameter :: neco=12, nsed=2  ! 
        integer(4),dimension(neco,0:n3d3nm0),target :: eco3nm0
        integer(4),dimension(neco,0:n3d3nm1),target :: eco3nm1
        integer(4),dimension(neco,0:n3d3nm2),target :: eco3nm2
        integer(4),dimension(neco,0:n3d3nm3),target :: eco3nm3
        integer(4),dimension(nsed,0:n2d3nm0),target :: sed3nm0
        integer(4),dimension(nsed,0:n2d3nm1),target :: sed3nm1
        integer(4),dimension(nsed,0:n2d3nm2),target :: sed3nm2
        integer(4),dimension(nsed,0:n2d3nm3),target :: sed3nm3
!       --------------------------------------------------
        real(4),allocatable :: buf_3d(:) 
        real(4),allocatable :: buf_2d(:)    
        integer(4),pointer  :: uu(:),vv(:),ww(:),ss(:),tt(:)  ! inherits size+offest
        integer(4),pointer  :: wu(:),wv(:),zz(:)              ! inherits size+offest
        integer(4),pointer  :: eco3nm(:,:)                    ! inherits size+offest
        integer(4),pointer  :: sed3nm(:,:)                    ! inherits size+offest
        character(len=2), parameter  :: area  = "ns"
        character*(*), parameter     :: usage = "Usage: burstHBM [p|b] <physfile|biofile> [prop*]"

        character cdate*10,ctmp*19,btmp*11
        character ncfname*999, prop*2, typ, action*24
        character*999 ::  ifilename
        real(4), parameter :: wuvscale = 1.0e7    ! wind scalefactor       -> m/s
        real(4), parameter :: uvscale  = 1.0e8    ! current scalefac.      -> m/s
        real(4), parameter :: stscale  = 1.0e7    ! salt/temp scfac.       -> PSU + deg.Celcius
        real(4), parameter :: wscale   = 1.0e11   ! vert.vel. scfac.       -> m
        real(4), parameter :: zscale   = 1.0e8    ! elevation scfac.       -> m
        real(4), parameter :: bscale   = 1.0e4    ! for biology conversion -> mmolN/m3

        integer :: ios, last, nfram, nargc, ia
        integer :: year,month,day,hour,minutes,seconds
        real(8) :: tmp
        logical :: readphys
!       ----------------------------------------------------------------------
!       select area
!       ----------------------------------------------------------------------
        select case (area)
           case ("ns")  ! -> area index 0
                uu     => uu0 
                vv     => vv0  
                ww     => ww0 
                ss     => ss0
                tt     => tt0
                wu     => wu0  
                wv     => wv0  
                zz     => zz0 
                eco3nm => eco3nm0
                sed3nm => sed3nm0
                allocate ( buf_3d(0:n3d3nm0) )  ! include element 0
                allocate ( buf_2d(0:n2d3nm0) )  ! include element 0
           case ("id") ! -> area index 1
                uu     => uu1 
                vv     => vv1  
                ww     => ww1 
                ss     => ss1
                tt     => tt1
                wu     => wu1  
                wv     => wv1  
                zz     => zz1 
                eco3nm => eco3nm1
                sed3nm => sed3nm1
                allocate ( buf_3d(0:n3d3nm1) )  ! include element 0
                allocate ( buf_2d(0:n2d3nm1) )  ! include element 0
           case ("ws") ! -> area index 2
                uu     => uu2 
                vv     => vv2  
                ww     => ww2 
                ss     => ss2
                tt     => tt2
                wu     => wu2  
                wv     => wv2  
                zz     => zz2 
                eco3nm => eco3nm2
                sed3nm => sed3nm2
                allocate ( buf_3d(0:n3d3nm2) )  ! include element 0
                allocate ( buf_2d(0:n2d3nm2) )  ! include element 0
           case ("bs") ! -> area index 3
                uu     => uu3 
                vv     => vv3  
                ww     => ww3 
                ss     => ss3
                tt     => tt3
                wu     => wu3  
                wv     => wv3  
                zz     => zz3 
                eco3nm => eco3nm3
                sed3nm => sed3nm3
                allocate ( buf_3d(0:n3d3nm3) )  ! include element 0
                allocate ( buf_2d(0:n2d3nm3) )  ! include element 0
        end select
!       ----------------------------------------------------------------------
!       Determine whether we read physics/biogeochemistry on first argument (p|b)
!       ----------------------------------------------------------------------
        nargc = iargc()
        if (nargc<2) then
           write(*,*) usage
           stop 2
        endif

        call getarg(1, typ)       ! 0 is program name
        if (typ == "p") then 
           readphys = .true.
           action   = "reading physics from : "
        elseif (typ == "b") then
           readphys = .false.
           action   = "reading biology from : "
        else
           write(*,*) "unrecognized file type:", typ
           stop 1
        endif

!       ----------------------------------------------------------------------
!       Read file name
!       ----------------------------------------------------------------------

        call getarg(2, ifilename)         ! 0 is program name
        open(11,file=ifilename,form='unformatted', status='old', iostat=ios)
        if(ios/=0) goto 800
        write(*,*) action, trim(adjustl(ifilename))
        
!       ----------------------------------------------------------------------
!       Loop over subsequent arguments (proterties to be extracted)
!       ----------------------------------------------------------------------   
     
        do ia = 3, nargc    ! possibly void loop
           call getarg(ia, prop) 
           write(*,*) "extracting ", prop, " from ", trim(adjustl(ifilename))
!          ------ read header
           rewind(11)
           read(11,end=850) btmp
           write(*,*) btmp
!          -------------- unconditional loop over frames (phys/bio specific) -----------
           nfram = 0
           do
              if (readphys) then
                 call read_cmod_frame(last)     ! last = 0: read single frame
              else
                 call read_ergom_frame(last)    ! last = 0: read single frame
              endif
              nfram = nfram + 1
           
              if (last==1) then    ! last = 1: EOF             
                 write(*,*) "no more frames in file"  
                 exit   ! frame loop, possibly continue property loop
              endif

              if (last==2) goto 850    ! last = 2: unexpected EOF
              
!         ------ dump frames depending on selection ----   
             
              select case (prop)

!             ------ physics ------
              case ("t")
                 if (.not.readphys) goto 891
                 write(ncfname, 534) "t_ns", year,month,day,hour,minutes,seconds
                 call get_t(buf_3d)
                 call dump_frame_as_netcdf_float4(trim(adjustl(ncfname)), buf_3d, 'deg.C')        

              case ("s")
                 if (.not.readphys) goto 891
                 write(ncfname, 534) "s_ns", year,month,day,hour,minutes,seconds
                 call get_s(buf_3d)
                 call dump_frame_as_netcdf_float4(trim(adjustl(ncfname)), buf_3d, 'PSU') 
                 
              case ("u")
                 if (.not.readphys) goto 891
                 write(ncfname, 534) "u_ns", year,month,day,hour,minutes,seconds
                 call get_u(buf_3d)
                 call dump_frame_as_netcdf_float4(trim(adjustl(ncfname)), buf_3d, 'east m/s') 

              case ("v")
                 if (.not.readphys) goto 891
                 write(ncfname, 534) "v_ns", year,month,day,hour,minutes,seconds
                 call get_v(buf_3d)
                 call dump_frame_as_netcdf_float4(trim(adjustl(ncfname)), buf_3d, 'north m/s') 

              case ("z")
                 if (.not.readphys) goto 891
                 write(ncfname, 534) "z_ns", year,month,day,hour,minutes,seconds
                 call get_z(buf_2d)   ! NB: 2D
                 call dump_frame_as_netcdf_float4(trim(adjustl(ncfname)), buf_2d, 'm') 

!             ------ biology ------
              case ("z1")
                 if (readphys) goto 890
                 write(ncfname, 534) "zo1_ns", year,month,day,hour,minutes,seconds
                 call get_zo1(buf_3d)
                 call dump_frame_as_netcdf_float4(trim(adjustl(ncfname)), buf_3d, 'mmol.m-3') 

              case ("z2")   
                 if (readphys) goto 890
                 write(ncfname, 534) "zo2_ns", year,month,day,hour,minutes,seconds
                 call get_zo2(buf_3d)
                 call dump_frame_as_netcdf_float4(trim(adjustl(ncfname)), buf_3d, 'mmol.m-3') 

              case ("ch")
                 if (readphys) goto 890
                 write(ncfname, 534) "chl_ns", year,month,day,hour,minutes,seconds
                 call get_chl(buf_3d)
                 call dump_frame_as_netcdf_float4(trim(adjustl(ncfname)), buf_3d, 'mmol.m-3') 

              case default
                 goto 880
                 
                 
           end select
           enddo  ! frame loop
        enddo     ! loop over arguments

        close(11)
        goto 900  
534     format(a,'_',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'_',i2.2,'_',i2.2,".nc")

!       ----------------------------------------------------------------------
!       write exit messages
!       ---------------------------------------------------------------------- 

800     write(*,*) "iostat error"
        goto 900
850     write(*,*) "unexpected EOF"
        goto 900
860     write(*,*) "illegal time stamp", ctmp
        goto 900
880     write(*,*) "unknown property request: ", prop
        goto 900
890     write(*,*) "inappropriate property request for physics: ", prop
        goto 900
891     write(*,*) "inappropriate property request for biology: ", prop
        goto 900

900     continue ! normal exit point



!       ===================================================================================
        contains
!       ===================================================================================
        


        subroutine read_cmod_frame(last) 
!       -------------------------------------------------------------
!       read single 4-area physics time frame from unit 11
!       uses main name space
!       last = 0: read single frame
!       last = 1: EOF
!       last = 2: unexpected EOF
!       -------------------------------------------------------------
        integer :: last
!       ------ time label ----     
        read(11,end=700) ctmp,tmp
        read(ctmp,532) year,month,day,hour,minutes,seconds   ! variables in global scope
532     format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)
        write(*,*) "read cmod frame ", ctmp, tmp
!       ------ data frames ----          
        read(11,end=800) wu0,wv0,wu1,wv1, wu2,wv2,wu3,wv3
        read(11,end=800) uu0,vv0,ww0,zz0
        read(11,end=800) uu1,vv1,ww1,zz1
        read(11,end=800) uu2,vv2,ww2,zz2
        read(11,end=800) uu3,vv3,ww3,zz3
        read(11,end=800) ss0,tt0
        read(11,end=800) ss1,tt1
        read(11,end=800) ss2,tt2
        read(11,end=800) ss3,tt3  
        last = 0  ! normal frame read
        return 
700     last = 1  ! EOF
        return 
800     last = 2  ! unexpected EOF
        return 
        end subroutine read_cmod_frame




        subroutine read_ergom_frame(last) 
!       -------------------------------------------------------------
!       read single 4-area ERGOM time frame from unit 11
!       uses main name space
!       last = 0: read single frame
!       last = 1: EOF
!       last = 2: unexpected EOF
!       -------------------------------------------------------------
        integer :: last
!       ------ time label ----     
        read(11,end=700) ctmp,tmp
        read(ctmp,532) year,month,day,hour,minutes,seconds
532     format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)
        write(*,*) "read ergom frame ", ctmp, tmp
!       ------ data frames ----          
        read(11,end=800) eco3nm0
        read(11,end=800) sed3nm0
        read(11,end=800) eco3nm1
        read(11,end=800) sed3nm1
        read(11,end=800) eco3nm2
        read(11,end=800) sed3nm2
        read(11,end=800) eco3nm3
        read(11,end=800) sed3nm3
        last = 0  ! normal frame read
        return 
700     last = 1  ! EOF
        return 
800     last = 2  ! unexpected EOF
        return 
        end subroutine read_ergom_frame


! ==============================================================================================
!       transform variables to float in main name space for current area 
!       result in buffer buf
! ==============================================================================================

!       ------ CMOD ------

        subroutine get_t(buf)
        real(4), intent(out) :: buf(:)
        buf = float(tt)/stscale    ! include element 0 in array, explicit cast to real*4
        end subroutine get_t

        subroutine get_s(buf)
        real(4), intent(out) :: buf(:)
        buf = float(ss)/stscale    ! include element 0 in array, explicit cast to real*4
        end subroutine get_s

        subroutine get_u(buf)
        real(4), intent(out) :: buf(:)
        buf = float(uu)/uvscale    ! include element 0 in array, explicit cast to real*4
        end subroutine get_u

        subroutine get_v(buf)
        real(4), intent(out) :: buf(:)
        buf = float(vv)/uvscale    ! include element 0 in array, explicit cast to real*4
        end subroutine get_v

        subroutine get_w(buf)
!       ------ transform t for current area ------
        real(4), intent(out) :: buf(:)
        buf = float(ww)/wscale    ! include element 0 in array, explicit cast to real*4
        end subroutine get_w

        subroutine get_z(buf)
        real(4), intent(out) :: buf(:)
        buf = float(zz)/zscale    ! include element 0 in array, explicit cast to real*4  
        end subroutine get_z

        subroutine get_wu(buf)
        real(4), intent(out) :: buf(:)
        buf = float(wu)/wuvscale    ! include element 0 in array, explicit cast to real*4  
        end subroutine get_wu

        subroutine get_wv(buf)
        real(4), intent(out) :: buf(:)
        buf = float(wv)/wuvscale    ! include element 0 in array, explicit cast to real*4  
        end subroutine get_wv

!       ------ ERGOM ------

        subroutine get_nh4(buf)
        real(4), intent(out) :: buf(:)
        buf = eco3nm(1,:)/bscale    ! include element 0 in array, explicit cast to real*4       
        end subroutine get_nh4
      
        subroutine get_no3(buf)
        real(4), intent(out) :: buf(:)
        buf = eco3nm(2, :)/bscale   ! include element 0 in array, explicit cast to real*4        
        end subroutine get_no3

        subroutine get_po4(buf)
        real(4), intent(out) :: buf(:)
        buf = eco3nm(3, :)/bscale   ! include element 0 in array, explicit cast to real*4         
        end subroutine get_po4
        
        subroutine get_Sio4(buf)
        real(4), intent(out) :: buf(:)
        buf = eco3nm(4, :)/bscale   ! include element 0 in array, explicit cast to real*4        
        end subroutine get_Sio4
             
        subroutine get_dia(buf)
        real(4), intent(out) :: buf(:)
        buf = eco3nm(5, :)/bscale  ! include element 0 in array, explicit cast to real*4         
        end subroutine get_dia

        subroutine get_fla(buf)
        real(4), intent(out) :: buf(:)
        buf = eco3nm(6, :)/bscale  ! include element 0 in array, explicit cast to real*4        
        end subroutine get_fla

        subroutine get_cya(buf)
        real(4), intent(out) :: buf(:)
        buf = eco3nm(7, :)/bscale  ! include element 0 in array, explicit cast to real*4         
        end subroutine get_cya

        subroutine get_zo1(buf)
        real(4), intent(out) :: buf(:)
        buf = eco3nm(8, :)/bscale   ! include element 0 in array, explicit cast to real*4        
        end subroutine get_zo1

        subroutine get_zo2(buf)
        real(4), intent(out) :: buf(:)
        buf = eco3nm(9, :)/bscale   ! include element 0 in array, explicit cast to real*4        
        end subroutine get_zo2

        subroutine get_SPM1(buf)
        real(4), intent(out) :: buf(:)
        buf = eco3nm(10, :)/bscale  ! include element 0 in array, explicit cast to real*4         
        end subroutine get_SPM1

        subroutine get_SPM2(buf)
        real(4), intent(out) :: buf(:)
        buf = eco3nm(11, :)/bscale  ! include element 0 in array, explicit cast to real*4         
        end subroutine get_SPM2

        subroutine get_oxy(buf)
        real(4), intent(out) :: buf(:)
        buf = eco3nm(12, :)/bscale  ! include element 0 in array, explicit cast to real*4         
        end subroutine get_oxy


        subroutine get_SPM_sed1(buf)
        real(4), intent(out) :: buf(:)
        buf =  sed3nm(1, :)/bscale  ! include element 0 in array, explicit cast to real*4         
        end subroutine get_SPM_sed1
      
        subroutine get_SPM_sed2(buf)
        real(4), intent(out) :: buf(:)
        buf =  sed3nm(2, :)/bscale   ! include element 0 in array, explicit cast to real*4        
        end subroutine get_SPM_sed2
       
!       ------ derived variables ------

        subroutine get_chl(buf)
        real(4), intent(out) :: buf(:)
        buf = (eco3nm(5,:) + eco3nm(6,:) + eco3nm(7,:))/bscale         ! include element 0 in array, explicit cast to real*4  
        end subroutine get_chl
        



        subroutine dump_frame_as_netcdf_float4(FILE_NAME, var, unit)
!       -------------------------------------------------------------
!       Dump real*4 array to netCDF file FILE_NAME
!       -------------------------------------------------------------
        character*(*), intent(in) :: FILE_NAME
        real(4), intent(in)       :: var(:)
        character*(*), intent(in) :: unit

        integer                   :: ncid, dimids(1)
        integer                   :: varid_data, varid_scale
!       -------------------------------------------------------------      
        call nccheck( nf90_create(FILE_NAME, NF90_CLOBBER, ncid)) 
 
!       --- create dims 
        call nccheck( nf90_def_dim(ncid, "nwet", size(var), dimids(1)) )  
            
!       --- create vars
        call nccheck( nf90_def_var(ncid, "data",  NF90_FLOAT,  dimids, varid_data) )  

!       --- create meta data
        call nccheck( nf90_put_att(ncid, NF90_GLOBAL, "unit", unit) )

!       --- end def mode 
        call nccheck( nf90_enddef(ncid) )

!       --- write variables
        call nccheck( nf90_put_var(ncid, varid_data, var) )

!       --- close set
        call nccheck( nf90_close(ncid) )

        end subroutine dump_frame_as_netcdf_float4


        subroutine nccheck(status)
!       -------------------------------------------------------------
!       Check status of netCDF operation and stop at error conditions
!       Apply as:
!           call nccheck( < netcdf_lib_call > )
!       -------------------------------------------------------------
        integer, intent(in) :: status
        if(status /= nf90_noerr) then 
           print *, trim(nf90_strerror(status))
           stop 3
        end if
        end subroutine nccheck  



        end program
