package ec.app.TracableProblems.TracableStatistics.ListOperations;

import ec.Individual;
import ec.vector.TracableDataTypes.TraceTuple;
import ec.vector.TracableDataTypes.TraceableString;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class IndividualAndGenomeListOperations {

    public static double getBest(ArrayList<Individual> inds){
        double highestFitness = -Double.MAX_VALUE;
        for(Individual i : inds)
            if (i.fitness.fitness() > highestFitness)
                highestFitness = i.fitness.fitness();
        return highestFitness;
    }
    public static double getMean(ArrayList<Individual> inds){
        double meanFitness = 0;
        for(Individual i : inds)
            meanFitness += i.fitness.fitness();
        meanFitness = meanFitness / inds.size();
        return meanFitness;
    }
    public static double getMedian(ArrayList<Individual> inds){
        ArrayList<Double> fitness = new ArrayList<Double>();
        for(Individual i : inds)
            fitness.add(i.fitness.fitness());
        Collections.sort(fitness);
        return fitness.get((int)(fitness.size() / 2));
    }


    /**
     * returns the genomes as a String array
     */
    public static String[] getGenotypeString(Individual ind){
        String[] array = ind.genotypeToString().split("\"");
        List<String> list = new ArrayList<String>(Arrays.asList(array));

        while (list.contains("")) {
            list.remove("");
        }

        list.remove(0); //remove the encoded genomeLength

        return list.toArray(new String[0]);
    }

    /**
     * returns the min fitness in a individual List
     */
    private static double getMinFitness(ArrayList<Individual> inds){
        double min = inds.get(0).fitness.fitness();

        for(int i = 1; i < inds.size(); i++)
            min = Math.min(inds.get(i).fitness.fitness(),min);

        return min;
    }

    /**
     * returns the max fitness in a individual List
     */
    private static double getMaxFitness(ArrayList<Individual> inds){
        double max = inds.get(0).fitness.fitness();

        for(int i = 1; i < inds.size(); i++)
            max = Math.max(inds.get(i).fitness.fitness(),max);

        return max;
    }

    /**
     * returns the genome length of Individual ind
     */
    public static int getGenomeLength(Individual ind){
        return Integer.parseInt(ind.genotypeToString().split("\"")[0].replaceAll("\\D+",""));
    }

    /**
     * Calculates the EntropyString, starting with the overall entropy of the generation following the entropy of every genome
     * @param genomes the genomes of the current generation as a TraceableString
     * @return the entropy string of the current generation as "overall entropy, entropy sum, entropy gene 1, entropy gene 2, ..."
     */
    public static String getEntropyString(ArrayList<ArrayList<TraceableString>> genomes){

        String entropyOfTheGeneration = "";

        int genomeCount = genomes.get(0).size();//getGenomeLength(inds.get(0));
        double entropySum = 0.0;

        for(int i = 0; i < genomeCount; i++){//iterate over all genomes (not individuals) //TODO: genomes.cout or genomeSize
            ArrayList<Integer> genomesToCheck = new ArrayList<Integer>();
            genomesToCheck.add(i);
            double currentEntropy = calculateEntropyOnGenomes(genomesToCheck, genomes);//TraceVectorListOperations.getShannonEntropy(traceVectors);
            entropyOfTheGeneration += ", " + currentEntropy;//calculateShannonEntropy(genomeEntropy);
            entropySum += currentEntropy;
        }

        ArrayList<Integer> genomesToCheck = new ArrayList<Integer>();
        for(int i = 0; i < genomeCount; i++) {//for to compute the overall entropy
            genomesToCheck.add(i);
        }
        return calculateEntropyOnGenomes(genomesToCheck, genomes) + ", " + entropySum + entropyOfTheGeneration;
    }


    //TODO: move to genomesListOperations
    /**
     * Returns the Entropy of a genome Column
     */
    private static double calculateEntropyOnGenomes(ArrayList<Integer> genomeColumns, ArrayList<ArrayList<TraceableString>> genomes){

        ArrayList<ArrayList<TraceTuple>> traceVectors = new ArrayList<ArrayList<TraceTuple>>(); //the traceVectors in genomeColumn
        for(int i = 0; i < genomeColumns.size(); i++) {
            int currentGenomeColumn = genomeColumns.get(i);
            if(currentGenomeColumn > genomes.get(0).size()) //NOTE: only works with fixed genome sizes
                throw new IndexOutOfBoundsException();
            for (int j = 0; j < genomes.size(); j++) {
                ArrayList<TraceableString> currentGenome = genomes.get(j);
                ArrayList<TraceTuple> currentTraceVector = currentGenome.get(currentGenomeColumn).getTraceVector();
                traceVectors.add(currentTraceVector);
            }
        }

        return TraceVectorListOperations.getShannonEntropy(traceVectors);
    }

    /**
     * returns a double Array with the length of three, containing in order:
     * countingImpact, fitnessImpact and entropyImpact
     * @param currentTraceID the current traceID to be checked
     * @param inds the individuals the traceID should be calculated for
     * @return [countingImpact, fitnessImpact, entropImpact, fitnessEntropyImpact]
     */
    public static double[] calculateImpact(int currentTraceID, ArrayList<Individual> inds, ArrayList<ArrayList<TraceableString>> genomes){ //TODO: delete k, input a entropy string?
        //TODO: input the TraceableString List instead of a ArrayList
        double countingImpact = 0;
        double fitnessImpact = 0;
        double entropyImpact = 0;
        double fitnessEntropyImpact = 0;
        double bestFitness = getMaxFitness(inds);
        double worstFitness = getMinFitness(inds);

        for(int i = 0; i < inds.size(); i++){//iterate over the individuals
            int popSize = inds.size();
            int genomeLength = getGenomeLength(inds.get(i));
            ArrayList<TraceableString> currentGenome = genomes.get(i);
            for(int j = 0; j < currentGenome.size(); j++){ //iterate over the genes
                TraceableString currentGene = currentGenome.get(j);
                List<TraceTuple> currentTraceVector = currentGene.getTraceVector();
                for(int k = 0; k < currentTraceVector.size(); k++){ //iterate over the traceVector
                    int traceID = currentTraceVector.get(k).getTraceID();
                    double influence = currentTraceVector.get(k).getImpact(); //the influence of the traceID in the gene
                    if(traceID == currentTraceID){
                        ArrayList<Integer> genomesToCheck = new ArrayList<>();
                        genomesToCheck.add(j);
                        double entropyFactor = calculateEntropyOnGenomes(genomesToCheck, genomes);
                        double diffToWorst = inds.get(i).fitness.fitness() - worstFitness;

                        countingImpact += influence/(popSize * genomeLength);
                        fitnessImpact += (influence/(popSize * genomeLength)) * (1 + diffToWorst);
                        entropyImpact += (influence/(popSize * genomeLength)) * (1 + entropyFactor);
                        fitnessEntropyImpact += (influence/(popSize * genomeLength)) * (1 + diffToWorst) * (1 + entropyFactor);
                    }
                }
            }
        }

        double[] returnArray = new double[4];
        returnArray[0] = countingImpact;
        returnArray[1] = fitnessImpact;
        returnArray[2] = entropyImpact;
        returnArray[3] = fitnessEntropyImpact;
        return returnArray;
    }

}
