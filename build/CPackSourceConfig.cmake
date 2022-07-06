# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. The list of available CPACK_xxx variables and their associated
# documentation may be obtained using
#  cpack --help-variable-list
#
# Some variables are common to all generators (e.g. CPACK_PACKAGE_NAME)
# and some are specific to a generator
# (e.g. CPACK_NSIS_EXTRA_INSTALL_COMMANDS). The generator specific variables
# usually begin with CPACK_<GENNAME>_xxxx.


set(CPACK_BUILD_SOURCE_DIRS "/Users/zhaoruiyang/libecc;/Users/zhaoruiyang/libecc/build")
set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
set(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
set(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
set(CPACK_DEBIAN_FILE_NAME "DEB-DEFAULT")
set(CPACK_DEBIAN_PACKAGE_NAME "libtommath-dev")
set(CPACK_DEBIAN_PACKAGE_RELEASE "1")
set(CPACK_DEBIAN_PACKAGE_SECTION "devel")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS "ON")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_FILE "/opt/homebrew/Cellar/cmake/3.22.3/share/cmake/Templates/CPack.GenericDescription.txt")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_SUMMARY "LIBECC built using CMake")
set(CPACK_FREEBSD_PACKAGE_CATEGORIES "math")
set(CPACK_FREEBSD_PACKAGE_MAINTAINER "gahr@FreeBSD.org")
set(CPACK_FREEBSD_PACKAGE_ORIGIN "math/libtommath")
set(CPACK_GENERATOR "TBZ2;TGZ;TXZ;TZ")
set(CPACK_IGNORE_FILES "/CVS/;/\\.svn/;/\\.bzr/;/\\.hg/;/\\.git/;\\.swp\$;\\.#;/#")
set(CPACK_INSTALLED_DIRECTORIES "/Users/zhaoruiyang/libecc;/")
set(CPACK_INSTALL_CMAKE_PROJECTS "")
set(CPACK_INSTALL_PREFIX "/usr/local")
set(CPACK_MODULE_PATH "")
set(CPACK_NSIS_DISPLAY_NAME "libtommath-devel 1.2.0")
set(CPACK_NSIS_INSTALLER_ICON_CODE "")
set(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
set(CPACK_NSIS_PACKAGE_NAME "libtommath-devel 1.2.0")
set(CPACK_NSIS_UNINSTALL_NAME "Uninstall")
set(CPACK_OSX_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX12.3.sdk")
set(CPACK_OUTPUT_CONFIG_FILE "/Users/zhaoruiyang/libecc/build/CPackConfig.cmake")
set(CPACK_PACKAGE_CONTACT "libtom@googlegroups.com")
set(CPACK_PACKAGE_DEFAULT_LOCATION "/")
set(CPACK_PACKAGE_DESCRIPTION_FILE "/opt/homebrew/Cellar/cmake/3.22.3/share/cmake/Templates/CPack.GenericDescription.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "LibTomMath")
set(CPACK_PACKAGE_DIRECTORY "/Users/zhaoruiyang/libecc/build/packages/Darwin/")
set(CPACK_PACKAGE_FILE_NAME "libtommath-devel-1.2.0-Source")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "libtommath-devel 1.2.0")
set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "libtommath-devel 1.2.0")
set(CPACK_PACKAGE_NAME "libtommath-devel")
set(CPACK_PACKAGE_RELOCATABLE "true")
set(CPACK_PACKAGE_VENDOR "libtom projects")
set(CPACK_PACKAGE_VERSION "1.2.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "2")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_RESOURCE_FILE_LICENSE "/Users/zhaoruiyang/libecc/libtommath/LICENSE")
set(CPACK_RESOURCE_FILE_README "/opt/homebrew/Cellar/cmake/3.22.3/share/cmake/Templates/CPack.GenericDescription.txt")
set(CPACK_RESOURCE_FILE_WELCOME "/opt/homebrew/Cellar/cmake/3.22.3/share/cmake/Templates/CPack.GenericWelcome.txt")
set(CPACK_RPM_PACKAGE_ARCHITECTURE "arm64")
set(CPACK_RPM_PACKAGE_LICENSE "The Unlicense")
set(CPACK_RPM_PACKAGE_NAME "libtommath-devel-1.2.0")
set(CPACK_RPM_PACKAGE_RELEASE "1")
set(CPACK_RPM_PACKAGE_SOURCES "ON")
set(CPACK_SET_DESTDIR "OFF")
set(CPACK_SOURCE_GENERATOR "TBZ2;TGZ;TXZ;TZ")
set(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;/\\.bzr/;/\\.hg/;/\\.git/;\\.swp\$;\\.#;/#")
set(CPACK_SOURCE_INSTALLED_DIRECTORIES "/Users/zhaoruiyang/libecc;/")
set(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/Users/zhaoruiyang/libecc/build/CPackSourceConfig.cmake")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "libtommath-devel-1.2.0-Source")
set(CPACK_SOURCE_RPM "OFF")
set(CPACK_SOURCE_TBZ2 "ON")
set(CPACK_SOURCE_TGZ "ON")
set(CPACK_SOURCE_TOPLEVEL_TAG "Darwin-Source")
set(CPACK_SOURCE_TXZ "ON")
set(CPACK_SOURCE_TZ "ON")
set(CPACK_SOURCE_ZIP "OFF")
set(CPACK_STRIP_FILES "")
set(CPACK_SYSTEM_NAME "Darwin")
set(CPACK_THREADS "1")
set(CPACK_TOPLEVEL_TAG "Darwin-Source")
set(CPACK_WIX_SIZEOF_VOID_P "8")

if(NOT CPACK_PROPERTIES_FILE)
  set(CPACK_PROPERTIES_FILE "/Users/zhaoruiyang/libecc/build/CPackProperties.cmake")
endif()

if(EXISTS ${CPACK_PROPERTIES_FILE})
  include(${CPACK_PROPERTIES_FILE})
endif()
