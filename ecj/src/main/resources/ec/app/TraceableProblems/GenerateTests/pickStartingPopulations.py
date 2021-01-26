import os
import pandas as pd

inputPath="../TestOutput/IntegerKnapsack/Hard/" #TODO: change
outputPath="../PickedPops/IntegerKnapsack/Hard/" #TODO: change

os.system("mkdir "+outputPath)



#read the dataframe
startingPopulationsDF = pd.read_csv(inputPath+"statisticOfGeneration.csv")


#choose the generations/copy the folder
one = startingPopulationsDF["bestFitness"].min()
five = startingPopulationsDF["bestFitness"].max()
two = one + ((five - one)/4)
three = one + ((five - one)/2)
four = one + (3*(five - one)/4)

picks = [startingPopulationsDF.iloc[(startingPopulationsDF['bestFitness']-one).abs().argsort()[:1]].index.values[0],
        startingPopulationsDF.iloc[(startingPopulationsDF['bestFitness']-two).abs().argsort()[:1]].index.values[0],
        startingPopulationsDF.iloc[(startingPopulationsDF['bestFitness']-three).abs().argsort()[:1]].index.values[0],
        startingPopulationsDF.iloc[(startingPopulationsDF['bestFitness']-four).abs().argsort()[:1]].index.values[0],
        startingPopulationsDF.iloc[(startingPopulationsDF['bestFitness']-five).abs().argsort()[:1]].index.values[0]]



count=0

for i in picks:
	currentInPath=inputPath+startingPopulationsDF.iloc[i]['folderName']+"/"
	currentOutPath=outputPath+"StartingPop"+str(count)+"/"
	os.system("mkdir "+currentOutPath)
	os.system("cp "+currentInPath+"startingPop.stat "+currentOutPath)
	count += 1
