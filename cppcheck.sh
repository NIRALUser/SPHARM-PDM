~/src/cppcheck/cppcheck \
    -i '/Users/johnsonhj/src/BRAINS4-build/BRAINSCommonLib-build/ITKIOFactoryRegistration /Users/johnsonhj/src/BRAINS4-build/include/ITK-4.0 /Users/johnsonhj/src/BRAINS4-extsrc/BRAINSCommonLib /Users/johnsonhj/src/BRAINS4-build/BRAINSCommonLib-build' \
--force --inline-suppr --template '{file}:{line},{severity},{id},{message}' \
--enable=all \
  $( find /Users/johnsonhj/src/spharm-pdm -name "*.h") \
                                      2> tee  err.txt

#--errorlist  \
#--xml  \

