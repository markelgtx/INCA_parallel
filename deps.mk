obj/inca_main.o : src/inca_main.f90 obj/mod_inputdat.o obj/mod_intrainfo.o obj/mod_input_reader.o obj/mod_read_files.o obj/mod_cubefile.o obj/mod_intracule.o obj/mod_c1hole.o
obj/mod_autogrid.o : src/mod_autogrid.f90 obj/mod_geninfo.o obj/mod_intrainfo.o
obj/mod_build_grid.o : src/mod_build_grid.f90 obj/mod_quadratures.o
obj/mod_c1hole.o : src/mod_c1hole.f90 obj/mod_geninfo.o obj/mod_wfxinfo.o obj/mod_numbers.o obj/mod_functions.o obj/mod_quadratures.o
obj/mod_constants.o : src/mod_constants.f90
obj/mod_cubefile.o : src/mod_cubefile.f90 obj/mod_geninfo.o obj/mod_wfxinfo.o obj/mod_functions.o obj/mod_numbers.o
obj/mod_functions.o : src/mod_functions.f90 obj/mod_geninfo.o obj/mod_wfxinfo.o obj/mod_loginfo.o obj/mod_numbers.o
obj/mod_geninfo.o : src/mod_geninfo.f90
obj/mod_input_reader.o : src/mod_input_reader.f90 obj/mod_inputdat.o obj/mod_intrainfo.o obj/mod_quadratures.o obj/mod_cubefile.o obj/mod_autogrid.o
obj/mod_inputdat.o : src/mod_inputdat.f90
obj/mod_intracule.o : src/mod_intracule.f90 obj/mod_geninfo.o obj/mod_numbers.o obj/mod_intrastuff.o obj/mod_wfxinfo.o obj/mod_intrainfo.o obj/mod_build_grid.o obj/mod_quadratures.o
obj/mod_intrainfo.o : src/mod_intrainfo.f90
obj/mod_intrastuff.o : src/mod_intrastuff.f90 obj/mod_geninfo.o obj/mod_quadratures.o obj/mod_numbers.o
obj/mod_locatemod.o : src/mod_locatemod.f90
obj/mod_loginfo.o : src/mod_loginfo.f90
obj/mod_numbers.o : src/mod_numbers.f90
obj/mod_quadratures.o : src/mod_quadratures.f90
obj/mod_radis.o : src/mod_radis.f90
obj/mod_read_files.o : src/mod_read_files.f90 obj/mod_wfxinfo.o obj/mod_geninfo.o obj/mod_loginfo.o obj/mod_locatemod.o obj/mod_numbers.o
obj/mod_wfxinfo.o : src/mod_wfxinfo.f90
