import csv
import os
import random


classpath = "../../../../../../../target/classes" #NOTE: change
executable = "ec.Evolve"

#statisticsFolder = "../TraceableVectorStatistics/TracableVectorStatisticsOut/"

#readCSVFile="generateEasyTests.csv"
print('path to the input file:')
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
def runOneTest(paramsFile, outputPath, numberRuns, additionalArguments, lineCount):
	os.makedirs(outputPath)
	appendCSVRow(outputPath+"statisticOfGeneration.csv", ["problem", "folderName", "seed"])
	for i in range(0, numberRuns):
		print("starting Run "+ str(i) +" of "+ str(numberRuns)+ " of problem " + str(lineCount))
		currentSeed=random.randint(-999999,999999)
		currentFolderName="Run"+str(i)
		currentRunPath=outputPath+currentFolderName+"/"
		os.makedirs(currentRunPath)
		evalFileParameter = " -p stat.eval-file="+ currentRunPath +"/eval.csv"
		startingPopFileParameter = " -p stat.starting-population-file="+ currentRunPath +"/startingPop.stat"
		seedParameter = " -p seed.0="+str(currentSeed)
		os.system("java -Xms1G -Xmx4G -cp "+classpath+" "+executable+" -file "+paramsFile + " "+ seedParameter + " " + evalFileParameter + " " + startingPopFileParameter + " " + additionalArguments)
		appendCSVRow(outputPath+"statisticOfGeneration.csv", [paramsFile, currentFolderName, str(currentSeed)])
		#os.system("cp -r "+statisticsFolder+" "+outputPath+currentFolderName )


#generate all the test in the csv file
csv_reader = csv.reader(open(readCSVFile,"r"),delimiter=',')
line_count = 0
for row in csv_reader:
	print("processing Problem ",line_count,"witn name", row[5])
	if line_count == 0:
		#headline of csv file
		print("Starting the Testruns!")
	else:
		if len(row) == 3:
			runOneTest(row[0], row[1], int(row[2]), "", line_count)
		if len(row) >= 4:
			runOneTest(row[0], row[1], int(row[2]), row[3], line_count)
	print("processed Problem "+str(line_count))
	line_count += 1
print(f'processed {line_count - 1} problems')
