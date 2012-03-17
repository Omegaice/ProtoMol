#ifndef PROTOMOL_PACKAGE_H
#define PROTOMOL_PACKAGE_H

#define PACKAGE_VERSION "3.3"
#define PACKAGE_BUGREPORT "protomol@cse.nd.edu, joseph@cauldrondevelopment.com"
#define PACKAGE_HOMEPAGE "http://protomol.sourceforge.net/"
#define PACKAGE_CITE \
  "T. Matthey, T. Cickovski, S. S. Hampton, A. Ko, Q. Ma,\n"\
  "M. Nyerges, T. Raeder, T. Slabach, and J. A. Izaguirre.\n"\
  "ProtoMol: An object-oriented framework for prototyping novel algorithms\n"\
  "for molecular dynamics.\n"\
  "ACM Trans. Math. Softw., 30(3):237 265, 2004."
#define PACKAGE_REVISION "${SVN_VERSION}"
#define PACKAGE_SOURCE_REPO "https://protomol.svn.sourceforge.net/svnroot/protomol/trunk/protomol"
#define PACKAGE_COMPILER "${CMAKE_CXX_COMPILER}"
#define PACKAGE_COMPILER_VERSION "4.2"
#define PACKAGE_COMPILER_FLAGS "${CMAKE_CXX_FLAGS}"
#define PACKAGE_COMPILER_LIBS "${LIBS}"
#define PACKAGE_BUILT_BY "${USER}"
#define PACKAGE_PLATFORM "${CMAKE_SYSTEM_NAME}"

#endif // PROTOMOL_PACKAGE_H
