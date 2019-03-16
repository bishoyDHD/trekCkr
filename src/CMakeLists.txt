CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)

add_subdirectory(tree)
add_subdirectory(InitReader)
add_subdirectory(RecipeReader)
add_subdirectory(plugins)
add_subdirectory(Plugin)
add_subdirectory(Chef)
add_subdirectory(cooker)
add_subdirectory(visualcooker)

add_subdirectory(SlowCtrl)

add_subdirectory(combine)



set (BUILD_WITHDB off CACHE bool "Enable DB link for database
(normally off)")

