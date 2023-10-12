Module ModuleFileName
    use mpi
    use M_strings_oop
    implicit none

    ! Global definition of file name extensions, paths, etc.
    integer(4), parameter :: FILE_EXTENSION_MODE_NONE = 11
    integer(4), parameter :: FILE_EXTENSION_MODE_TXT = 12
    integer(4), parameter :: FILE_EXTENSION_MODE_DAT = 13
    integer(4), parameter :: FILE_EXTENSION_MODE_JSON = 14
    integer(4), parameter :: FILE_EXTENSION_MODE_H5 = 15
    integer(4), parameter :: FILE_EXTENSION_MODE_MILO = 16

    integer(4), parameter :: FILE_PATH_MODE_NONE = 21
    integer(4), parameter :: FILE_PATH_MODE_RESTART = 22
    integer(4), parameter :: FILE_PATH_MODE_DIAG = 23
    integer(4), parameter :: FILE_PATH_MODE_CHECK = 24
    integer(4), parameter :: FILE_PATH_MODE_INPUT = 25
    integer(4), parameter :: FILE_PATH_MODE_INIT = 26

    integer(4), parameter :: FILE_PARALLEL_MODE_CLOSE = 31          ! 是否需要添加image的序号
    integer(4), parameter :: FILE_PARALLEL_MODE_OPEN = 32

    integer(4), parameter :: FILE_DYNAMIC_MODE_CLOSE = -1           ! 动态命名文件
    integer(4), parameter :: FILE_DYNAMIC_MODE_OPEN = 0


    type FileName
        type(string)    :: DataName                                 ! 文件标识名
        type(string)    :: FullName                                 ! 文件完整的相对路径

        integer(4)      :: PathMode         = FILE_PATH_MODE_NONE
        integer(4)      :: ExtensionMode    = FILE_EXTENSION_MODE_NONE
        integer(4)      :: ParallelMode     = FILE_PARALLEL_MODE_CLOSE
        integer(4)      :: DynamicIndex     = FILE_DYNAMIC_MODE_CLOSE

    contains
        procedure :: InitFromPara => InitializationFileName
        procedure :: InitFromtype => InitializationFileNameFromOthertype
        generic   :: Init => InitFromPara, InitFromtype

        procedure :: SetPath => SetFileNamePathMode
        procedure :: SetExte => SetFileNameExtensionMode
        procedure :: SetParl => SetFileNameParallelMode
        procedure :: SetIndex=>SetFileNameDynamicIndex

        procedure :: UpdateIndex => UpdateFileNameDynamicIndex

        procedure, private :: Update => UpdateFileNameFullName
   
    end type FileName


    ! Global group definition
    type(FileName), save :: RESTART_FILE_NAME
    type(FileName), save :: DIAG_FILE_NAME
    type(FileName), save :: CHECK_FILE_NAME


    contains
    
        subroutine InitFileName()

            RESTART_FILE_NAME%DataName%str = "restart"
            RESTART_FILE_NAME%PathMode = FILE_PATH_MODE_RESTART
            RESTART_FILE_NAME%ExtensionMode = FILE_EXTENSION_MODE_H5
            RESTART_FILE_NAME%ParallelMode = FILE_PARALLEL_MODE_OPEN
            RESTART_FILE_NAME%DynamicIndex = FILE_DYNAMIC_MODE_CLOSE

            DIAG_FILE_NAME%DataName%str = "diag"
            DIAG_FILE_NAME%PathMode = FILE_PATH_MODE_DIAG
            DIAG_FILE_NAME%ExtensionMode = FILE_EXTENSION_MODE_DAT
            DIAG_FILE_NAME%ParallelMode = FILE_PARALLEL_MODE_CLOSE
            DIAG_FILE_NAME%DynamicIndex = FILE_DYNAMIC_MODE_CLOSE

            CHECK_FILE_NAME%DataName%str = "check"
            CHECK_FILE_NAME%PathMode = FILE_PATH_MODE_CHECK
            CHECK_FILE_NAME%ExtensionMode = FILE_EXTENSION_MODE_DAT
            CHECK_FILE_NAME%ParallelMode = FILE_PARALLEL_MODE_CLOSE
            CHECK_FILE_NAME%DynamicIndex = FILE_DYNAMIC_MODE_CLOSE

        end subroutine InitFileName


        subroutine InitializationFileName(IOName, DataNameInp, PathMode, ExtensionMode, ParallelMode, DynamicIndex)
            class(FileName), intent(inout)  :: IOName
            character(*), intent(in) :: DataNameInp
            integer(4), intent(in), Optional :: PathMode, ExtensionMode, ParallelMode, DynamicIndex

            IOName%DataName%str = Trim(DataNameInp)

            if(present(PathMode))       call IOName%SetPath(PathMode)
            if(present(ExtensionMode))  call IOName%SetExte(ExtensionMode)
            if(present(ParallelMode))   call IOName%SetParl(ParallelMode)
            if(present(DynamicIndex))   call IOName%SetIndex(DynamicIndex)

            call IOName%Update()

        end subroutine InitializationFileName


        subroutine InitializationFileNameFromOthertype(IOName, DataNameInp, FileNameConst)
            class(FileName), intent(inout)  :: IOName
            character(*), intent(in) :: DataNameInp
            class(FileName), intent(in)  :: FileNameConst

            IOName%PathMode = FileNameConst%PathMode
            IOName%ExtensionMode = FileNameConst%ExtensionMode
            IOName%ParallelMode = FileNameConst%ParallelMode
            IOName%DynamicIndex = FileNameConst%DynamicIndex

            IOName%DataName%str = Trim(DataNameInp)

            call IOName%Update()

        end subroutine InitializationFileNameFromOthertype


        subroutine SetFileNamePathMode(IOName, PathMode)
            class(FileName), intent(inout)  :: IOName
            integer(4), intent(in) :: PathMode

            IOName%PathMode = PathMode
            call IOName%Update()
            
        end subroutine SetFileNamePathMode


        subroutine SetFileNameExtensionMode(IOName, ExtensionMode)
            class(FileName), intent(inout)  :: IOName
            integer(4), intent(in) :: ExtensionMode

            IOName%ExtensionMode = ExtensionMode
            call IOName%Update()

        end subroutine SetFileNameExtensionMode


        subroutine SetFileNameParallelMode(IOName, ParallelMode)
            class(FileName), intent(inout)  :: IOName
            integer(4), intent(in) :: ParallelMode

            IOName%ParallelMode = ParallelMode
            call IOName%Update()
            
        end subroutine SetFileNameParallelMode
        

        subroutine SetFileNameDynamicIndex(IOName, DynamicIndex)
            class(FileName), intent(inout)  :: IOName
            integer(4), intent(in) :: DynamicIndex

            IOName%DynamicIndex = DynamicIndex
            if (IOName%DynamicIndex < 0) IOName%DynamicIndex = FILE_DYNAMIC_MODE_CLOSE
            
            call IOName%Update()
            
        end subroutine SetFileNameDynamicIndex


        subroutine UpdateFileNameDynamicIndex(IOName)
            class(FileName), intent(inout)  :: IOName

            if (IOName%DynamicIndex < 0) then
                IOName%DynamicIndex = FILE_DYNAMIC_MODE_CLOSE
            else
                IOName%DynamicIndex = IOName%DynamicIndex + 1
            end if

            call IOName%Update()

        end subroutine UpdateFileNameDynamicIndex


        subroutine UpdateFileNameFullName(IOName)
            class(FileName),intent(inout)  :: IOName
            integer(4) :: size, rank, ierr

            call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

            associate(DataName      => IOName%DataName,     FullName        => IOName%FullName, &
                      PathMode      => IOName%PathMode,     ExtensionMode   => IOName%ExtensionMode, &
                      ParallelMode  => IOName%ParallelMode, DynamicIndex    => IOName%DynamicIndex)

                select case(PathMode)
                    case(FILE_PATH_MODE_NONE)
                        FullName%str = './'

                    case(FILE_PATH_MODE_RESTART)
                        FullName%str = './restart/'

                    case(FILE_PATH_MODE_DIAG)
                        FullName%str = './diag/'

                    case(FILE_PATH_MODE_CHECK)
                        FullName%str = './check/'

                    case(FILE_PATH_MODE_INPUT)
                        FullName%str = './input/'

                    case(FILE_PATH_MODE_INIT)
                        FullName%str = './init/'

                    case default
                        FullName%str = './'
                end select

                FullName%str = FullName%str // DataName%str

                select case(ParallelMode)
                    case(FILE_PARALLEL_MODE_OPEN)
                        FullName%str = FullName%str // '_' // trim(num2str(rank, 3))

                    case(FILE_PARALLEL_MODE_CLOSE)
                    case default
                end select

                select case(DynamicIndex)
                    case(FILE_DYNAMIC_MODE_CLOSE)

                    case(FILE_DYNAMIC_MODE_OPEN)
                        FullName%str = FullName%str // '_' // trim(num2str(DynamicIndex, 5))

                    case default
                        FullName%str = FullName%str // '_' // trim(num2str(DynamicIndex, 5))
                end select

                select case(ExtensionMode)
                    case(FILE_EXTENSION_MODE_TXT)
                        FullName%str = FullName%str // '.txt'

                    case(FILE_EXTENSION_MODE_DAT)
                        FullName%str = FullName%str // '.dat'

                    case(FILE_EXTENSION_MODE_H5)
                        FullName%str = FullName%str // '.h5'

                    case(FILE_EXTENSION_MODE_JSON)
                        FullName%str = FullName%str // '.json'

                    case(FILE_EXTENSION_MODE_MILO)
                        FullName%str = FullName%str // '.milo'

                    case(FILE_EXTENSION_MODE_NONE)
                    case default
                end select

            end associate

        end subroutine UpdateFileNameFullName


        function num2str(digit, width)
            character(len=99) :: num2str
            integer(4), intent(in) :: digit, width
            integer(4) :: digit_w, i, digit_tmp

            num2str = ''
            digit_w = getDigitWidth(digit)
            if (digit >= 0 .and. digit_w <= width) then

                digit_tmp = digit
                do i = width, 1, -1
                    if (digit_tmp > 0) then
                        num2str(i:i) = achar(mod(digit_tmp, 10) + iachar('0'))
                        digit_tmp = floor(dble(digit_tmp) / 10.d0)
                    else
                        num2str(i:i) = '0'
                    end if
                end do

            else
                do i = 1, width
                    num2str(i:i) = '0'
                end do
            end if
        end


        function getDigitWidth(digit)
            integer(4) :: getDigitWidth
            integer(4), intent(in) :: digit

            if (digit < 0) then
                getDigitWidth = 0

            else if (digit < 10) then
                getDigitWidth = 1

            else
                getDigitWidth = floor(log10(dble(digit))) + 1

            end if
        end

end Module ModuleFileName

