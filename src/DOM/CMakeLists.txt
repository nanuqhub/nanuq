set(MY_MODS daymod.F90
            domain.F90
            dombth.F90
            domhgr.F90
            dommsk.F90
            dom_oce.F90
            domutl.F90
            domwri.F90
            istate.F90
            phycst.F90
            CACHE INTERNAL "")
add_library(mylib "${MY_MODS}")
target_link_libraries(PRIVATE mylib)
