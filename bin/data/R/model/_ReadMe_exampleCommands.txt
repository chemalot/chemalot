
cd /gne/research/data/smdd/ApplicationData/ModelData/RModelTest/

#---Example 1: Create and predict RandomForest Models---#
set timestamp=`date +%Y%m%d%H%M`

  ~smdi/dev/bin/sdfRRandomForestCreator.pl -in testData/HLM_set_3k_TRAIN.sdf \
      -modelLocation RModelTest/manle/HLM_set_3k_$timestamp/RFModel \
      -modelName test_RF_train \
      -responseField measured_CL -descriptorFields "cLogP|MW" \
      -predictionColName "cHLM_CL" -identifierField G

  ~smdi/dev/bin/sdfRModelPredictor.pl -in testData/HLM_set_3k_TEST_diverse.sdf \
      -modelLocation RModelTest/manle/HLM_set_3k_$timestamp/RFModel \
      -modelName test_RF_train \
      -NNFieldPrefix cHLM \
      -out manle/HLM_set_3k_$timestamp/out/V1.out.sdf


#---Example 2: Create and predict RandomForest Models---#
set timestamp=`date +%Y%m%d%H%M`

  ~smdi/dev/bin/sdfRRandomForestCreator.pl -in testData/HLM_set_3k_TRAIN.sdf \
      -modelLocation RModelTest/manle/HLM_set_3k_$timestamp/RFModel \
      -modelName test_RF_train \
      -responseField measured_CL -descriptorFields "cLogP|MW" \
      -predictionColName "cHLM_CL" \
      -nTree 5 -mTry 5 -logImportance \
      -conversionScript testData/HLM_set_3k_conversionFunctions.R

  ~smdi/dev/bin/sdfRModelPredictor.pl -in testData/HLM_set_3k_TRAIN.sdf \
      -modelLocation RModelTest/manle/HLM_set_3k_$timestamp/RFModel \
      -modelName test_RF_train \
      -NNFieldPrefix cHLM \
      -out manle/HLM_set_3k_$timestamp/out/V1.out.sdf \
      -computeOutputScript testData/HLM_set_3k_computeOutput.R


#---Example 3: Create and predict SVM Models---#
set timestamp=`date +%Y%m%d%H%M`

  ~smdi/dev/bin/sdfRSVMCreator.pl -in testData/Test_zeroVariance_collinearity.tab \
      -modelLocation ./manle/HLM_set_3k_$timestamp/SVMModel  \
      -modelName test_collinear \
      -responseField measured_CL -descriptorFields "cLogP|MW|test_collinear|zero_variance" \
      -predictionColName "cHLM_CL" \
      -gamma "2/nx" -cost "1.0" -crossValidationFolds 5 -type Default \
      -kernel Radial \
      -conversionScript testData/HLM_set_3k_conversionFunctions.R
  
  ~smdi/dev/bin/sdfRModelPredictor.pl -in testData/HLM_set_3k_TRAIN.sdf \
      -modelLocation ./manle/HLM_set_3k_$timestamp/SVMModel \
      -modelName test_collinear \
      -NNSimAbove 0.5 -NNFieldPrefix cHLM \
      -out ./manle/HLM_set_3k_$timestamp/out/V1.out.sdf \
      -computeOutputScript testData/HLM_set_3k_computeOutput.R
# Error Message (Should be more clear that descriptor fields are missing in the 
#                test data set)
#   Expected to get as many predictions as inputs. inputs= 100 predictions= 0


#---Example 4: Create and predict SVM Models---#
set timestamp=`date +%Y%m%d%H%M`

  ~smdi/dev/bin/sdfRSVMCreator.pl -in testData/HLM_set_3k_TRAIN.sdf \
      -modelLocation ./manle/HLM_set_3k_$timestamp/SVMModel  \
      -modelName test_SVM_train \
      -responseField measured_CL -descriptorFields "cLogP|MW" \
      -predictionColName "cHLM_CL" \
      -identifierField G \
      -gamma "1/nx,2/nx" -cost "1.0,2.0" -crossValidationFolds 5 -type Default \
      -kernel Radial \
      -conversionScript testData/HLM_set_3k_conversionFunctions.R

  ~smdi/dev/bin/sdfRModelPredictor.pl -in testData/HLM_set_3k_TRAIN.sdf \
      -modelLocation ./manle/HLM_set_3k_$timestamp/SVMModel \
      -modelName test_SVM_train \
      -NNSimAbove 0.5 -NNFieldPrefix cHLM \
      -out ./manle/HLM_set_3k_$timestamp/out/V1.out.sdf \
      -computeOutputScript testData/HLM_set_3k_computeOutput.R

  ~smdi/dev/bin/sdfRModelPredictor.pl -in testData/HLM_set_3k_TEST_diverse.sdf \
      -modelLocation ./manle/HLM_set_3k_$timestamp/SVMModel \
      -modelName test_SVM_train \
      -NNSimAbove 0.545 -NNFieldPrefix cHLM \
      -out ./manle/HLM_set_3k_$timestamp/out/V2.out.sdf \
      -computeOutputScript testData/HLM_set_3k_computeOutput.R


#---Example 5: Test the removeNAandINF() function in sdfRCreatorFunctions.R---#
set timestamp=`date +%Y%m%d%H%M`

  ~smdi/dev/bin/sdfRRandomForestCreator.pl -in testData/ContainsNA.sdf \
      -modelLocation RModelTest/manle/ContainsNA_$timestamp/RFModel \
      -modelName ContainsNA \
      -responseField measured_CL -descriptorFields "cLogP|MW" \
      -predictionColName "cHLM_CL" -identifierField G

  ~smdi/dev/bin/sdfRSVMCreator.pl -in testData/ContainsNA.sdf \
      -modelLocation ./manle/ContainsNA_$timestamp/SVMModel  \
      -modelName ContainsNA \
      -responseField measured_CL -descriptorFields "cLogP|MW" \
      -predictionColName "cHLM_CL" \
      -gamma "2/nx" -cost "1.0" -crossValidationFolds 5 -type Default \
      -kernel Radial



#---Create RandomForest Models---#
  ~smdi/dev/bin/sdfRRandomForestCreator.pl -in testData/Test_zeroVariance_collinearity.tab \
      -modelLocation ./manle/HLM_set_3k_RFModel_3 -modelName test_collinear \
      -responseField measured_CL -descriptorFields "cLogP|MW|test_collinear|zero_variance" \
      -predictionColName "cHLM_CL" \
      -nTree 1 -mTry 1 -logImportance \
      -conversionScript testData/HLM_set_3k_conversionFunctions.R


#---Create SVM Models---#
  ~smdi/dev/bin/sdfRSVMCreator.pl -in testData/HLM_set_3k_TRAIN.sdf \
      -modelLocation RModelTest/manle/HLM_set_3k_SVMModel_1  \
      -modelName test_SVM_train \
      -responseField measured_CL -descriptorFields "cLogP|MW" -predictionColName "cHLM_CL" \
      -gamma "1/nx,2/nx" -cost "1.0,2.0" -crossValidationFolds 5 -type Default \
      -kernel Radial

  ~smdi/dev/bin/sdfRSVMCreator.pl -in testData/HLM_set_3k_TRAIN.sdf \
      -modelLocation RModelTest/manle/HLM_set_3k_SVMModel_2  \
      -modelName test_SVM_train \
      -responseField measured_CL -descriptorFields "cLogP|MW" \
      -predictionColName "cHLM_CL" \
      -gamma "1/nx,2/nx" -cost "1.0,2.0" -crossValidationFolds 5 -type Default \
      -kernel Radial \
      -conversionScript testData/HLM_set_3k_badConversionFunctions.R

  ~smdi/dev/bin/sdfRSVMCreator.pl -in testData/HLM_set_3k_TRAIN.sdf \
      -modelLocation RModelTest/manle/HLM_set_3k_SVMModel_2.only  \
      -modelName test_SVM_train \
      -responseField measured_CL -descriptorFields "cLogP|MW" \
      -predictionColName "cHLM_CL" \
      -gamma "1/nx,2/nx" -cost "1.0,2.0" -crossValidationFolds 5 -type Default \
      -kernel Radial \
      -conversionScript testData/HLM_set_3k_onlyConversionFunctions.R

  ~smdi/dev/bin/sdfRSVMCreator.pl -in testData/HLM_set_3k_TRAIN.sdf \
      -modelLocation RModelTest/manle/HLM_set_3k_SVMModel_2  \
      -modelName test_SVM_train \
      -responseField measured_CL -descriptorFields "cLogP|MW" \
      -predictionColName "cHLM_CL" \
      -gamma "1/nx,2/nx" -cost "1.0,2.0" -crossValidationFolds 5 -type Default \
      -kernel Radial \
      -conversionScript testData/HLM_set_3k_conversionFunctions.R





