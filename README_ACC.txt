Notes on WW3 ACC Directives 01/05/2024 
======================================
1. To get nvtx to link using cmake
(This worked on Glados at DL, there may be other ways):

Add the line:

target_link_libraries(ww3_lib PUBLIC nvhpcwrapnvtx)

in the file

model/src/CMakeLists.txt

2. The ACC directives in file w3srcemd.F90 are currently commented out.
When optimal CPU=>GPU data transfer points in the source code have been determined (for unmanaged memory) then these may be re-inserted.
