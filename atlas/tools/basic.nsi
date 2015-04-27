;===========================================================================
; This script will build an installer with only mingw32 .a static libraries
; See superpack.nsi for a complete installer with import libraries and dll.
;===========================================================================

!define PRODUCT_NAME "BlasLapack superpack"

;--------------------------------
;Extra headers necessary

!include 'MUI.nsh'
!include 'Sections.nsh'
!include 'LogicLib.nsh'

SetCompress off ; Useful to disable compression under development

;--------------------------------
;General

;Name and file
Name "${PRODUCT_NAME}"
OutFile "setup.exe"

;Default installation folder
InstallDir "$PROGRAMFILES\${PRODUCT_NAME}"

;Get installation folder from registry if available
InstallDirRegKey HKCU "Software\${PRODUCT_NAME}" ""

;--------------------------------
;Interface Settings

!define MUI_ABORTWARNING
#!addplugindir "/home/david/local/share/nsis/plugins"

;--------------------------------
;Pages

;!insertmacro MUI_PAGE_LICENSE "${NSISDIR}\Docs\Modern UI\License.txt"
!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES

!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
;  !insertmacro MUI_PAGE_FINISH_TEXT "Installation done ! DO NOT FORGET TO PUT THE INSTALLED DLL SOMEWHERE WHERE THEY CAN BE FOUND BY WINDOWS !"

;--------------------------------
;Languages

!insertmacro MUI_LANGUAGE "English"

;--------------------------
; Types of installation
; - default: install most optimized available, mingw only (toward people who just want thing to work)
; - custom: choose whatever you want
; - install everything
InstType "Typical"
InstType "Full"

;--------------------------------
;Component Sections

Var HasSSE2
Var HasSSE3
Var CPUSSE

SubSection "Generic BLAS" BLAS
        Section "Mingw archives" BLAS_MINGW
                SetOutPath "$INSTDIR\Generic\mingw32"
                File "binaries\nosse\mingw32\libblas.a"
        SectionEnd
SubSectionEnd

SubSection "Generic LAPACK" LAPACK
        Section "Mingw archives" LAPACK_MINGW
                SetOutPath "$INSTDIR\Generic\mingw32"
                File "binaries\nosse\mingw32\liblapack.a"
        SectionEnd
SubSectionEnd

SubSection "ATLAS" ATLAS
        SubSection "SSE2" ATLAS_SSE2
                Section "Mingw archives" ATLAS_SSE2_MINGW
                        SetOutPath "$INSTDIR\ATLAS\sse2\mingw32"
                        File "binaries\sse2\mingw32\lib*.a"
                SectionEnd
        SubSectionEnd

        SubSection "SSE3" ATLAS_SSE3
                Section "Mingw archives" ATLAS_SSE3_MINGW
                        SetOutPath "$INSTDIR\ATLAS\sse3\mingw32"
                        File "binaries\sse3\mingw32\lib*.a"
                SectionEnd
        SubSectionEnd
SubSectionEnd

Section -Post
        ;Create uninstaller
        WriteUninstaller "$INSTDIR\Uninstall.exe"
SectionEnd

Section -Dummy DUMMY
SectionEnd

Function .onInit
        ; Detect CPU target capabilities for ATLAS
        StrCpy $CPUSSE "0"
        CpuCaps::hasSSE3
        Pop $0
        StrCpy $HasSSE3 $0
        CpuCaps::hasSSE2
        Pop $0
        StrCpy $HasSSE2 $0

        ; Take care about the order !

        ; SSE 2
        StrCmp $HasSSE2 "Y" include_sse2 no_include_sse2
        include_sse2:
                DetailPrint '"Target CPU handles SSE2"'
                StrCpy $CPUSSE "2"
                goto done_sse2
        no_include_sse2:
                DetailPrint '"Target CPU doe NOT handle SSE2"'
                goto done_sse2
        done_sse2:

        ; SSE 3
        StrCmp $HasSSE3 "Y" include_sse3 no_include_sse3
        include_sse3:
                DetailPrint '"Target CPU handles SSE3"'
                StrCpy $CPUSSE "3"
                goto done_sse3
        no_include_sse3:
                DetailPrint '"Target CPU doe NOT handle SSE3"'
                goto done_sse3
        done_sse3:

        ${Switch} $CPUSSE
        ${Case} "3"
                DetailPrint '"Install SSE 3"'
                !insertmacro SetSectionInInstType   "${ATLAS_SSE3_MINGW}" "1"
                ${Break}
        ${Case} "2"
                DetailPrint '"Install SSE 2"'
                !insertmacro SetSectionInInstType   "${ATLAS_SSE2_MINGW}" "1"
                ${Break}
        ${Default}
                DetailPrint '"Install NO SSE"'
                !insertmacro SetSectionInInstType   "${BLAS_MINGW}" "1"
                !insertmacro SetSectionInInstType   "${LAPACK_MINGW}" "1"
                ${Break}
        ${EndSwitch}

        !insertmacro SetSectionInInstType "${BLAS_MINGW}" "2"

        !insertmacro SetSectionInInstType "${LAPACK_MINGW}" "2"

        !insertmacro SetSectionInInstType "${ATLAS_SSE2_MINGW}" "2"

        !insertmacro SetSectionInInstType "${ATLAS_SSE3_MINGW}" "2"

        SetCurInstType "0"
        WriteRegStr HKCU "Software\${PRODUCT_NAME}" "" $INSTDIR
FunctionEnd

;--------------------------------
;Descriptions
;Language strings
LangString DESC_BLAS ${LANG_ENGLISH} "This section contains BLAS libraries.\
        Those are generic, non optimized libraries (from netlib)."
LangString DESC_LAPACK ${LANG_ENGLISH} "This section contains LAPACK libraries.\
        Those are generic, non optimized libraries (from netlib)."

LangString DESC_ATLAS ${LANG_ENGLISH} "This section contains ATLAS libraries.\
        Those are optimized libraries for BLAS and LAPACK."
LangString DESC_ATLAS_SSE2 ${LANG_ENGLISH} "This section contains ATLAS libraries \
        for CPU supporting at least SSE2. "
LangString DESC_ATLAS_SSE3 ${LANG_ENGLISH} "This section contains ATLAS libraries \
        for CPU supporting at least SSE3. "

LangString DESC_MINGW ${LANG_ENGLISH} "This will install .a libraries: \
        this is what you need if you compile numpy with mingw32 compilers."

;Assign language strings to sections
!insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
        !insertmacro MUI_DESCRIPTION_TEXT ${BLAS} $(DESC_BLAS)
        !insertmacro MUI_DESCRIPTION_TEXT ${BLAS_MINGW} $(DESC_MINGW)

        !insertmacro MUI_DESCRIPTION_TEXT ${LAPACK} $(DESC_LAPACK)
        !insertmacro MUI_DESCRIPTION_TEXT ${LAPACK_MINGW} $(DESC_MINGW)

        !insertmacro MUI_DESCRIPTION_TEXT ${ATLAS} $(DESC_ATLAS)
        !insertmacro MUI_DESCRIPTION_TEXT ${ATLAS_SSE2} $(DESC_ATLAS_SSE2)
        !insertmacro MUI_DESCRIPTION_TEXT ${ATLAS_SSE2_MINGW} $(DESC_MINGW)

        !insertmacro MUI_DESCRIPTION_TEXT ${ATLAS_SSE3} $(DESC_ATLAS_SSE3)
        !insertmacro MUI_DESCRIPTION_TEXT ${ATLAS_SSE3_MINGW} $(DESC_MINGW)
!insertmacro MUI_FUNCTION_DESCRIPTION_END

;--------------------------------
;Uninstaller Section

Section "Uninstall"

        Delete "$INSTDIR\Uninstall.exe"

        ; XXX: we should track the installed files instead, but I am lazy

        ; Generic libs
        Delete "$INSTDIR\generic\mingw32\libblas.a"
        Delete "$INSTDIR\generic\mingw32\liblapack.a"
        RMDir "$INSTDIR\generic\mingw32"

        RMDir "$INSTDIR\generic\lib"
        RMDir "$INSTDIR\generic"

        ; ATLAS SSE2
        Delete "$INSTDIR\sse2\mingw32\*.a"
        RMDir "$INSTDIR\sse2\mingw32"
        RMDir "$INSTDIR\sse2"

        ; ATLAS SSE2
        Delete "$INSTDIR\sse3\mingw32\*.a"
        RMDir "$INSTDIR\sse3\mingw32"
        RMDir "$INSTDIR\sse3"
        RMDir "$INSTDIR"

        DeleteRegKey /ifempty HKCU "Software\${PRODUCT_NAME}"

SectionEnd
