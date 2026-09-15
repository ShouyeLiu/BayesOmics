add_test(NAME public_cli_help COMMAND $<TARGET_FILE:BayesOmics64> --help)
add_test(NAME public_cli_reject_sampler
  COMMAND ${CMAKE_COMMAND} -DBINARY=$<TARGET_FILE:BayesOmics64>
    -DCASE=sampler -P ${CMAKE_CURRENT_SOURCE_DIR}/check-public-cli.cmake)
add_test(NAME public_cli_reject_model
  COMMAND ${CMAKE_COMMAND} -DBINARY=$<TARGET_FILE:BayesOmics64>
    -DCASE=model -P ${CMAKE_CURRENT_SOURCE_DIR}/check-public-cli.cmake)
