# Short overview of the python scripts and how to use them

The scripts "generateFromCSV.py" and "pickStartingpopulations.py" help with running more than one test configuration at once. 

## generateFromCSV.py

There is a little setup-process with the script when using it on another computer. You need to change a few Paths inside of the script:
* classpath: your Java-Classpath
* executable: path to the ECJ-Executalbe (ec.Evolve)
* statisticsFolder: path to the statistic-output of the ECJ-Problem(s) you wand to run (relative od absolute)

Having set the correct Paths, the script can be run. Running the Script, you are promted to input a (path to a) CSV-File, including the informations aboud the tests that should be run. The CSV-File should be structured like this:

 * the first row is just used for headlining the columns, they will not be read from the script.
 * first column: first program-argument (-file ...), relative or absolute path to the .params file
 * second column: The outputPath for the test-result. Make shure that this path does not exist jet! (otherwise the scipt will stop with an error)
 * third column: The number of runs per test configuration
 * fourth column: optional column for additional arguments, like the Seed (-seed ...) or the StartingPopulation (-p ...)
 * the following columns are optional and contain information that could be used in visualization scripts:
 	* fifth column: specifies the problem type for the visualization
 	* sixth column: specifies the difficulty for the visualization scripts with Easy, Medium and Hard
 	* seventh column: specifies the name of the population for the visualization

Example.csv shows an example configuration file.
