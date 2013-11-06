# for documentation please look at files in ./Pro2_doc/
#
# for a free and powerful f90 compiler please check out gfortran 
# Prospino2.1 is now tested on Mac OS-X 10.5 gfortran (gcc-4.3.0)
#
# f90 compilers which are tested: Intel, Dec, Portland Group, NAG, gfortran
#
# switch off optimization for Intel compiler v.8, otherwise VEGAS will die
#OPTION = -O0
# these options needed for the NAG compiler
#OPTION = -ieee=full -gline -dcfuns#
#
# paths, compiler, etc -> feel free to edit 
#
COMP   = gfortran
DIRECT = ./
#
# no editing should be required below, unless you want to change the pdf set
#
MACROS = ${DIRECT}Pro2_subroutines/
INTEGS = ${DIRECT}Pro2_integrals/
MATRIX = ${DIRECT}Pro2_matrix/
INTERF = ${DIRECT}Pro2_interface/
STRONG = ${DIRECT}Pro2_sq_gl/
#
# sets of files
#
FILES_DIRECT = ${DIRECT}Xvital.o \
               ${DIRECT}Xin_out.o \
               ${DIRECT}Xinteg_ng.o \
               ${DIRECT}Xinteg_ns.o \
               ${DIRECT}Xinteg_nn.o \
               ${DIRECT}Xinteg_ll.o \
               ${DIRECT}Xinteg_tb.o \
               ${DIRECT}Xinteg_lq.o \
               ${DIRECT}Xinteg_le.o \
               ${DIRECT}Xinteg_hh.o \
               ${DIRECT}Xinteg_ht.o \
               ${DIRECT}Xinitialize.o \
               ${DIRECT}Xprospino_subroutine.o
FILES_MACROS = ${MACROS}Xsugra.o \
               ${MACROS}Xalphas_smooth.o \
               ${MACROS}Xrunning_mass.o \
               ${MACROS}Xvegas_array.o \
	       ${MACROS}Xread_les_houches.o \
	       ${MACROS}Cteq6Pdf-2008.o 
#tq5	       ${MACROS}XCtq5Pdf.o 
#tq6	       ${MACROS}XCteq61Pdf.o 
#tqx	       ${MACROS}Cteq6Pdf-2008.o 
FILES_INTEGS = ${INTEGS}Xtwopoint.o \
	       ${INTEGS}Xthreepoint.o \
	       ${INTEGS}Xfourpoint.o \
	       ${INTEGS}Xangular_basic.o \
	       ${INTEGS}Xangular_functions_lq.o \
	       ${MACROS}Xangular_array_ng.o \
	       ${MACROS}Xangular_array_ns.o \
	       ${MACROS}Xangular_array_tb.o \
	       ${MACROS}Xangular_array_lq.o \
	       ${MACROS}Xangular_array_le.o \
	       ${MACROS}Xangular_array_hh.o \
	       ${MACROS}Xangular_array_ht.o \
	       ${MACROS}Xsoft_array_lq.o \
	       ${MACROS}Xscalar_array_ng.o \
	       ${MACROS}Xscalar_array_ns.o \
	       ${MACROS}Xscalar_array_nn.o \
	       ${MACROS}Xscalar_array_tb.o \
	       ${MACROS}Xscalar_array_lq.o \
	       ${MACROS}Xscalar_array_le.o \
	       ${MACROS}Xscalar_array_hh.o \
	       ${MACROS}Xscalar_array_ht.o 
FILES_MATRIX = ${MATRIX}Xmatrix_ng_v.o \
	       ${MATRIX}Xmatrix_ng_r.o \
	       ${MATRIX}Xmatrix_ng_c.o \
	       ${MATRIX}Xmatrix_ns_v.o \
	       ${MATRIX}Xmatrix_ns_qg.o \
	       ${MATRIX}Xmatrix_ns_gg.o \
	       ${MATRIX}Xmatrix_ns_qb.o \
	       ${MATRIX}Xmatrix_ns_qq.o \
	       ${MATRIX}Xmatrix_nn_v.o \
	       ${MATRIX}Xmatrix_nn_r.o \
	       ${MATRIX}Xmatrix_nn_c.o \
	       ${MATRIX}Xmatrix_ll_v.o \
	       ${MATRIX}Xmatrix_ll_r.o \
	       ${MATRIX}Xmatrix_ll_c.o \
	       ${MATRIX}Xmatrix_tb_qb.o \
	       ${MATRIX}Xmatrix_tb_gg.o \
	       ${MATRIX}Xmatrix_tb_cr.o \
	       ${MATRIX}Xmatrix_bb_qb.o \
	       ${MATRIX}Xmatrix_bb_gg.o \
	       ${MATRIX}Xmatrix_lq_qb.o \
	       ${MATRIX}Xmatrix_lq_gg.o \
	       ${MATRIX}Xmatrix_lq_cr.o \
	       ${MATRIX}Xmatrix_le_v.o \
	       ${MATRIX}Xmatrix_le_qg.o \
	       ${MATRIX}Xmatrix_le_gg.o \
	       ${MATRIX}Xmatrix_le_qb.o \
	       ${MATRIX}Xmatrix_le_qq.o \
	       ${MATRIX}Xmatrix_hh_v.o \
	       ${MATRIX}Xmatrix_hh_gb.o \
	       ${MATRIX}Xmatrix_hh_qb.o \
	       ${MATRIX}Xmatrix_hh_qg.o \
	       ${MATRIX}Xmatrix_hh_s.o \
	       ${MATRIX}Xmatrix_ht_v.o \
	       ${MATRIX}Xmatrix_ht_s.o \
	       ${MATRIX}Xmatrix_ht_qq.o \
	       ${MATRIX}Xmatrix_ht_qb.o \
	       ${MATRIX}Xmatrix_ht_qg.o \
	       ${MATRIX}Xmatrix_ht_gg.o \
	       ${MACROS}Xlumis.o
FILES_INTERF = ${INTERF}Xget_pdf.o \
	       ${INTERF}Xget_spectrum.o \
	       ${INTERF}Xhard_stop.o
FILES_STRONG = ${STRONG}hadrongg.o \
	       ${STRONG}hadronsb.o \
	       ${STRONG}hadronsg.o \
	       ${STRONG}hadronss.o \
	       ${STRONG}hadronst.o \
	       ${STRONG}initpdf.o \
	       ${STRONG}integral.o \
	       ${STRONG}integralst.o \
	       ${STRONG}matrixgg.o \
	       ${STRONG}matrixsb.o \
	       ${STRONG}matrixsg.o \
	       ${STRONG}matrixss.o \
	       ${STRONG}matrixst.o 
#
# build the prospino executable
#
prospino:	prospino_main.f90  ${FILES_DIRECT} \
	                           ${FILES_MACROS} \
                                   ${FILES_MATRIX} \
                                   ${FILES_INTEGS} \
                                   ${FILES_STRONG} \
                                   ${FILES_INTERF} 
	${COMP} prospino_main.f90  ${FILES_DIRECT} \
	                           ${FILES_MACROS} \
                                   ${FILES_MATRIX} \
                                   ${FILES_INTEGS} \
                                   ${FILES_STRONG} \
                                   ${FILES_INTERF} -o prospino_2.run
#
#
#
validate:	     validate.f90  ${FILES_DIRECT} \
	                           ${FILES_MACROS} \
                                   ${FILES_MATRIX} \
                                   ${FILES_INTEGS} \
                                   ${FILES_STRONG} \
                                   ${FILES_INTERF} 
	${COMP}      validate.f90  ${FILES_DIRECT} \
	                           ${FILES_MACROS} \
                                   ${FILES_MATRIX} \
                                   ${FILES_INTEGS} \
                                   ${FILES_STRONG} \
                                   ${FILES_INTERF} -o validate.run
#
#
#
clean:
	rm -i ${DIRECT}*.mod ; \
	rm -i $(DIRECT)prospino_2.run ; \
	rm -i $(DIRECT)validate.run ; \
	rm -i ${DIRECT}*.o ; rm -i ${DIRECT}*~ ; \
	rm -i ${MACROS}*.o ; rm -i ${MACROS}*~ ; \
	rm -i ${INTEGS}*.o ; rm -i ${INTEGS}*~ ; \
	rm -i ${MATRIX}*.o ; rm -i ${MATRIX}*~ ; \
	rm -i ${STRONG}*.o ; rm -i ${STRONG}*~ ; \
	rm -i ${INTERF}*.o ; rm -i ${INTERF}*~   
#
# prepare all the object files from the main directory [DIRECT]
# for DEC compiler keep them in order of dependence 
#
${DIRECT}%.o:  ${DIRECT}%.f90
	${COMP} -c -o $@ $<
#
# prepare all the object files from the directory 'Pro2_interface' [INTERF]
#
${INTERF}%.o:  ${INTERF}%.f
	${COMP} -c -o $@ $<
#
# prepare all the object files from the directory 'Pro2_subroutines' [MACROS]
#
${MACROS}%.o:  ${MACROS}%.f
	${COMP} -c -o $@ $<
#
# prepare all the object files from the directory 'Pro2_integrals' [INTEGS]
#
${INTEGS}%.o:  ${INTEGS}%.f
	${COMP} -c -o $@ $<
#
# prepare all the object files from the directory 'Pro2_matrix' [MATRIX]
#
${MATRIX}%.o:  ${MATRIX}%.f
	${COMP} -c -o $@ $<
#
# prepare all the object files from the directory 'Pro2_sq_gl' [STRONG]
#
${STRONG}%.o:  ${STRONG}%.f
	${COMP} -c -o $@ $<
#



