# Included by project() of a subproject that should not get link-time optimization:
# the builtin LLVM, Clang and Cling of ROOT. They are not in the hot path, and LTO
# of them costs much build time and memory.
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
