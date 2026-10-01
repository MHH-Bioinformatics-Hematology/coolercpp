# The programs the documentation shows. They are built with the library so that
# an example that stops compiling fails the build rather than the reader.
foreach(example read_cool create_cool scool_fileops)
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE coolercpp::coolercpp)
endforeach()
