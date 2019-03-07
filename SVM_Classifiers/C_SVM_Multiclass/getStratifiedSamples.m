function [trnDataIndxArr, testDataIndxArr] = getStratifiedSamples(trainDataperClass, dataPerClass, noOfClasses)
    trnDataIndxArr=[];
    testDataIndxArr=[];
    for i = 1:noOfClasses
        startIndx = (dataPerClass*(i-1))+1;
        endIndex = dataPerClass*i;
        trainDataForClass = randperm(dataPerClass,trainDataperClass) + (dataPerClass*(i-1)) ;
        testDataForClass = setdiff(startIndx:endIndex,trainDataForClass);
        trnDataIndxArr = [trnDataIndxArr trainDataForClass];
        testDataIndxArr = [testDataIndxArr testDataForClass];
    end      
end