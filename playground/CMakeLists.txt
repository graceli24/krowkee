
function(setup_krowkee_target exe_name)
	add_executable(${exe_name} ${exe_name}.cpp)
    target_link_libraries(${exe_name} PRIVATE krowkee)
	message(STATUS "Adding include directory: " ${CMAKE_SOURCE_DIR}/include)
	target_include_directories(${exe_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
endfunction()

setup_krowkee_target(int128_test)
setup_krowkee_target(kpoly_hash_test)