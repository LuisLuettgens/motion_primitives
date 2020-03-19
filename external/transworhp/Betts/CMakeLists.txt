project(bettsTestSet CXX Fortran)

set(PROBS_PATH problems)

add_library(${PROJECT_NAME}_problems SHARED

	auxiliaryFunctions.f
	${PROBS_PATH}/bdprob.f

	${PROBS_PATH}/alprde.f		${PROBS_PATH}/alprin.f	${PROBS_PATH}/spyke.f
	${PROBS_PATH}/aquade.f		${PROBS_PATH}/aquain.f
	${PROBS_PATH}/araode.f		${PROBS_PATH}/araoin.f
	${PROBS_PATH}/ashrde.f		${PROBS_PATH}/ashrin.f
#	${PROBS_PATH}/bangde.f		${PROBS_PATH}/bangin.f	${PROBS_PATH}/linkph.f
	${PROBS_PATH}/bracde.f		${PROBS_PATH}/bracin.f
	${PROBS_PATH}/brgrde.f		${PROBS_PATH}/brgrin.f	${PROBS_PATH}/brgreq.f
	${PROBS_PATH}/chmrde.f		${PROBS_PATH}/chmrin.f	${PROBS_PATH}/chmrdf.f
	${PROBS_PATH}/crande.f		${PROBS_PATH}/cranin.f
	${PROBS_PATH}/ffrbde.f		${PROBS_PATH}/ffrbin.f
#	${PROBS_PATH}/gydnde.f		${PROBS_PATH}/gydnin.f	${PROBS_PATH}/gydndf.f ${PROBS_PATH}/gydngs.f fehlt: atmml, hdenvl, hdmvgi
	${PROBS_PATH}/hangde.f		${PROBS_PATH}/hangin.f # freie Endzeit.. aber nicht .f ???
	${PROBS_PATH}/heatde.f		${PROBS_PATH}/heatin.f	${PROBS_PATH}/pderhs.f ${PROBS_PATH}/qsubi.f
	${PROBS_PATH}/jshide.f		${PROBS_PATH}/jshiin.f
	${PROBS_PATH}/kplrde.f		${PROBS_PATH}/kplrin.f
	${PROBS_PATH}/lnhtde.f		${PROBS_PATH}/lnhtin.f
	${PROBS_PATH}/lntsde.f		${PROBS_PATH}/lntsin.f # Randbedinung =45 fehlt..
	${PROBS_PATH}/lowtde.f		${PROBS_PATH}/lowtin.f
	${PROBS_PATH}/medide.f		${PROBS_PATH}/mediin.f
	${PROBS_PATH}/pndlde.f		${PROBS_PATH}/pndlin.f
	${PROBS_PATH}/qlinde.f		${PROBS_PATH}/qlinin.f # funktion fehlt "hdmvps_"
#	${PROBS_PATH}/raylde.f		${PROBS_PATH}/raylin.f funktionen fehlen..
	${PROBS_PATH}/rbrmde.f		${PROBS_PATH}/rbrmin.f
#	${PROBS_PATH}/rivrde.f		${PROBS_PATH}/rivrin.f fehlt: hdfrvl, pthcon, autopl
#	${PROBS_PATH}/robode.f		${PROBS_PATH}/roboin.f fehlt: usrdeq, dircom, linkst
#	${PROBS_PATH}/soarde.f		${PROBS_PATH}/soarin.f fehlt: pthcon, autolk
	${PROBS_PATH}/tb2sde.f		${PROBS_PATH}/tb2sin.f
	${PROBS_PATH}/trande.f		${PROBS_PATH}/tranin.f
	${PROBS_PATH}/zrmlde.f		${PROBS_PATH}/zrmlin.f
)

target_link_libraries(${PROJECT_NAME}_problems
	PRIVATE
		blas
)

add_executable (${PROJECT_NAME} testSet.cpp)

target_include_directories(
	${PROJECT_NAME}
	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}"
		$<TARGET_PROPERTY:TransWORHP_core,INTERFACE_INCLUDE_DIRECTORIES>
)

target_link_libraries(
	${PROJECT_NAME}
	PUBLIC
		TransWORHP
		${PROJECT_NAME}_problems
)

target_compile_options(
	${PROJECT_NAME}
	PUBLIC
		-std=c++11
	PRIVATE
		-Wall
		-Wextra
		-pedantic
)

