import csv
import os
import random


classpath = "../../../../../../../target/classes" #NOTE: change
executable = "ec.Evolve"
#statisticsFolder="../TracableStatistics/JustStartingPopulationStatisticsOut/" 
statisticsFolder="../TracableStatistics/TracableVectorStatisticsOut/" #NOTE: change

#readCSVFile="generateEasyTests.csv"
readCSVFile=input()

if(len(readCSVFile) == 0):
	exit()

def readCSVRow (filename, rowNumber):
	csv_reader = csv.reader(open(filename,"r"),delimiter=',')
	rowCount=0
	for row in csv_reader:
		if rowCount == rowNumber:
			return row
		else:
			rowCount += 1


def appendCSVRow (fileName, csvRow):
    csv_writer = csv.writer(open(fileName, mode='a'), delimiter=',')
    csv_writer.writerow(csvRow)


#define a method for running one test
def runOneTest(paramsFile, outputPath, numberRuns, additionalArguments):
	os.makedirs(outputPath)
	appendCSVRow(outputPath+"statisticOfGeneration.csv", ["problem", "folderName", "seed", "bestFitness", "meanFitness", "medianFitness"])
	for i in range(0, numberRuns):
		print("starting Run "+ str(i) +" of "+ str(numberRuns))
		seed=random.randint(-999999,999999)
		currentFolderName="Run"+str(i)
		os.system("java -cp "+classpath+" "+executable+" -file "+paramsFile+" -p seed.0="+str(seed) + " " + additionalArguments)
		startingEndingFitness=readCSVRow(statisticsFolder+"startingEndingFitness.csv", 1)
		appendCSVRow(outputPath+"statisticOfGeneration.csv", [paramsFile, currentFolderName, str(seed), startingEndingFitness[0], startingEndingFitness[1], startingEndingFitness[2]])		
		os.system("cp -r "+statisticsFolder+" "+outputPath+currentFolderName )






#generate all the test in the csv file
csv_reader = csv.reader(open(readCSVFile,"r"),delimiter=',')
line_count = 0
for row in csv_reader:
	print("processing Problem "+str(line_count))
	if line_count == 0:
		#headline of csv file
		print("Starting the Testruns!")
	else:
		if len(row) == 3:
			runOneTest(row[0], row[1], int(row[2]), "")
		if len(row) >= 4:
			runOneTest(row[0], row[1], int(row[2]), row[3])
	print("processed Problem "+str(line_count))
	line_count += 1
print(f'processed {line_count - 1} problems')
