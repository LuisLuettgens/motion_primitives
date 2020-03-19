file(REMOVE_RECURSE
  "libTransWORHP.pdb"
  "libTransWORHP.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/TransWORHP.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
