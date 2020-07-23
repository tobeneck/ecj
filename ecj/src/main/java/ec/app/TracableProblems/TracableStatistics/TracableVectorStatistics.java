package ec.app.TracableProblems.TracableStatistics;

import ec.EvolutionState;
import ec.Individual;
import ec.Statistics;
import ec.app.TracableProblems.TracableStatistics.ListOperations.DoubleListOperations;
import ec.app.TracableProblems.TracableStatistics.ListOperations.IndividualListOperations;
import ec.util.Parameter;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class TracableVectorStatistics extends Statistics
{
    // The parameter string and log number of the file for our readable population
    public static final String P_POPFILE = "pop-file";
    public int popLog;

    public static final String P_STATFILE = "stat-file"; //TODO: rename
    public int statLog;

    public static final String P_STARTINGPOP = "starting-population-file";
    public int startingPopulationLog;

    public static final String P_FITNESSFILE = "fitness-file";
    public int fitnessLog;

    public static final String P_TRACEIDSFILE= "traceIDs-file";
    public int traceIDsLog;

    public static final String P_GENOMEFILE= "genome-file";
    public int genomeLog;

    public static final String P_COUNTINGIMPACTFILE = "counting-impact-file";
    public int countingImpactLog;

    public static final String P_FITNESSIMPACTFILE = "fitness-impact-file";
    public int fitnessImpactLog;

    public static final String P_ENTROPYIMPACTFILE = "entropy-impact-file";
    public int entropyImpactLog;

    public static final String P_FITNESSENTROPYIMPACTFILE = "fitness-entropy-impact-file";
    public int fitnessEntropyImpactLog;

    public static final String P_NORMALIZEDFITNESSIMPACTFILE = "normalized-fitness-impact-file";
    public int normalizedFitnessImpactLog;

    public static final String P_NORMALIZEDENTROPYIMPACTFILE = "normalized-entropy-impact-file";
    public int normalizedEntropyImpactLog;

    public static final String P_NORMALIZEDFITNESSENTROPYIMPACTFILE = "normalized-fitness-entropy-impact-file";
    public int normalizedFitnessEntropyImpactLog;

    public static final String P_ENTROPYFILE = "entropy-file";
    public int entropyLog;



    String fitnessString = "";

    String traceIDsString = "";
    String genomeString = "";

    String impactHeadder = "";
    String countingImpactString = "";
    String fitnessImpactString = "";
    String entropyImpactString = "";
    String fitnessEntropyImpactString = "";
    String normalizedFitnessImpactString = "";
    String normalizedEntropyImpactString = "";
    String normalizedFitnessEntropyImpactString = "";
    String entropyString = "";
    String startingEndingFitness = "";


    public void setup(final EvolutionState state, final Parameter base)
    {
        // DO NOT FORGET to call super.setup(...) !!
        super.setup(state,base);

        // set up popFile
        File popFile = state.parameters.getFile(
                base.push(P_POPFILE),null);
        if (popFile!=null) try
        {
            popLog = state.output.addLog(popFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    popFile + ":\n" + i);
        }


        // set up statFile
        File statFile = state.parameters.getFile(
                base.push(P_STATFILE),null);
        if (statFile!=null) try
        {
            statLog = state.output.addLog(statFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    statFile + ":\n" + i);
        }

        // set up startingPopulationFile
        File startingPopulationFile = state.parameters.getFile(
                base.push(P_STARTINGPOP),null);
        if (startingPopulationFile!=null) try
        {
            startingPopulationLog = state.output.addLog(startingPopulationFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    startingPopulationFile + ":\n" + i);
        }

        // set up fitnessFile
        File fitnessFile = state.parameters.getFile(
                base.push(P_FITNESSFILE),null);
        if (fitnessFile!=null) try
        {
            fitnessLog = state.output.addLog(fitnessFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    fitnessFile + ":\n" + i);
        }

        // set up traceIDsFile
        File traceIDsFile = state.parameters.getFile(
                base.push(P_TRACEIDSFILE),null);
        if (traceIDsFile!=null) try
        {
            traceIDsLog = state.output.addLog(traceIDsFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    traceIDsFile + ":\n" + i);
        }

        // set up genomeFile
        File genomeFile = state.parameters.getFile(
                base.push(P_GENOMEFILE),null);
        if (genomeFile!=null) try
        {
            genomeLog = state.output.addLog(genomeFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    genomeFile + ":\n" + i);
        }

        // set up CountinIgmpactFile
        File countingImpactFile = state.parameters.getFile(
                base.push(P_COUNTINGIMPACTFILE),null);
        if (countingImpactFile!=null) try
        {
            countingImpactLog = state.output.addLog(countingImpactFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    countingImpactFile + ":\n" + i);
        }

        // set up fitnessImpactFile
        File fitnessImpactFile = state.parameters.getFile(
                base.push(P_FITNESSIMPACTFILE),null);
        if (fitnessImpactFile!=null) try
        {
            fitnessImpactLog = state.output.addLog(fitnessImpactFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    fitnessImpactFile + ":\n" + i);
        }

        // set up entropyImpactFile
        File entropyImpactFile = state.parameters.getFile(
                base.push(P_ENTROPYIMPACTFILE),null);
        if (entropyImpactFile!=null) try
        {
            entropyImpactLog = state.output.addLog(entropyImpactFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    entropyImpactFile + ":\n" + i);
        }

        // set up fitnessEntropyImpactFile
        File fitnessEntropyImpactFile = state.parameters.getFile(
                base.push(P_FITNESSENTROPYIMPACTFILE),null);
        if (fitnessEntropyImpactFile!=null) try
        {
            fitnessEntropyImpactLog = state.output.addLog(fitnessEntropyImpactFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    fitnessEntropyImpactFile + ":\n" + i);
        }

        // set up normalizedFitnessImpactFile
        File normalizedFitnessImpactFile = state.parameters.getFile(
                base.push(P_NORMALIZEDFITNESSIMPACTFILE),null);
        if (normalizedFitnessImpactFile!=null) try
        {
            normalizedFitnessImpactLog = state.output.addLog(normalizedFitnessImpactFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    normalizedFitnessImpactFile + ":\n" + i);
        }

        // set up normalizedEntropyImpactFile
        File normalizedEntropyImpactFile = state.parameters.getFile(
                base.push(P_NORMALIZEDENTROPYIMPACTFILE),null);
        if (normalizedEntropyImpactFile!=null) try
        {
            normalizedEntropyImpactLog = state.output.addLog(normalizedEntropyImpactFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    normalizedEntropyImpactFile + ":\n" + i);
        }

        // set up normalizedFitnessEntropyImpactFile
        File normalizedFitnessEntropyImpactFile = state.parameters.getFile(
                base.push(P_NORMALIZEDFITNESSENTROPYIMPACTFILE),null);
        if (normalizedFitnessEntropyImpactFile!=null) try
        {
            normalizedFitnessEntropyImpactLog = state.output.addLog(normalizedFitnessEntropyImpactFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    normalizedFitnessEntropyImpactFile + ":\n" + i);
        }

        // set up entropyFile
        File entropyFile = state.parameters.getFile(
                base.push(P_ENTROPYFILE),null);
        if (entropyFile!=null) try
        {
            entropyLog = state.output.addLog(entropyFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    entropyFile + ":\n" + i);
        }

    }



    public void postEvaluationStatistics(final EvolutionState state)
    {
        // be certain to call the hook on super!
        //super.postEvaluationStatistics(state);

        ArrayList<Individual> inds = state.population.subpops.get(0).individuals;

        //print the startingPopulationFile and the startingFitness
        if(state.generation == 0){
            // print out the population
            state.population.subpops.get(0).printSubpopulation(state,startingPopulationLog);
            startingEndingFitness = "" + IndividualListOperations.getBest(inds) + "," + IndividualListOperations.getMean(inds) +","+IndividualListOperations.getMedian(inds);
        }


        //print the pop.stat
        // write out a warning that the next generation is coming
        state.output.println("-----------------------\nGENERATION " +
                state.generation + "\n-----------------------", popLog);
        // print out the population
        state.population.printPopulation(state,popLog);


        //print the impactFiles
        String fitnessOut = state.generation +"";
        for(int i = 0; i < inds.size(); i++) {
            fitnessOut += "," + inds.get(i).fitness.fitness();
        }
        fitnessString += fitnessOut + "\n";

        //print the traceIDsFile
        for(int i = 0; i < inds.size(); i++) {
            String out = state.generation + "," + i + ",";
            String[] genotypeString = inds.get(i).genotypeToString().split("\"");
            for(int j = 1; j < genotypeString.length; j++){
                if(!genotypeString[j].isEmpty())
                    out = out + genotypeString[j].split(",")[1] + ",";
            }

            out = out.substring(0, out.length() - 1);
            traceIDsString += out + "\n";
        }

        //print the genomeFile
        for(int i = 0; i < inds.size(); i++) {
            String out = state.generation + "," + i + ",";
            String[] genotypeString = inds.get(i).genotypeToString().split("\"");
            for(int j = 1; j < genotypeString.length; j++){
                if(!genotypeString[j].isEmpty())
                    out = out + genotypeString[j].split(",")[0] + ",";
            }

            out = out.substring(0, out.length() - 1);
            genomeString += out + "\n";
        }

        //go over all initial genomes
        List<Double> countingImpact = new ArrayList<Double>();
        List<Double> fitnessImpact = new ArrayList<Double>();
        List<Double> entropyImpact = new ArrayList<Double>();
        List<Double> fitnessEntropyImpact = new ArrayList<Double>();
        for(int i = 0; i < inds.size(); i++) {
            double[] impact = IndividualListOperations.calculateImpact(i, inds);
            countingImpact.add(impact[0]);
            fitnessImpact.add(impact[1]);
            entropyImpact.add(impact[2]);
            fitnessEntropyImpact.add(impact[3]);
        }

        //go over all mutations again for mutation, seperate to build a sum
        for(int i = -1; i >= state.mutationCounter; i--) {
            double[] impact = IndividualListOperations.calculateImpact(i, inds);
            countingImpact.add(impact[0]);
            fitnessImpact.add(impact[1]);
            entropyImpact.add(impact[2]);
            fitnessEntropyImpact.add(impact[3]);
        }

        //print the impactFile
        List<Double> countingImpactHead = countingImpact.subList(0,inds.size());
        List<Double> fitnessImpactHead = fitnessImpact.subList(0,inds.size());
        List<Double> entropyImpactHead = entropyImpact.subList(0,inds.size());
        List<Double> fitnessEntropyImpactHead = entropyImpact.subList(0,inds.size());
        List<Double> countingImpactTail = new ArrayList<Double>();
        List<Double> fitnessImpactTail = new ArrayList<Double>();
        List<Double> entropyImpactTail = new ArrayList<Double>();
        List<Double> fitnessEntropyImpactTail = new ArrayList<Double>();
        if(state.mutationCounter < 0){
            countingImpactTail = countingImpact.subList(inds.size(),countingImpact.size());
            fitnessImpactTail = fitnessImpact.subList(inds.size(),fitnessImpact.size());
            entropyImpactTail = entropyImpact.subList(inds.size(),entropyImpact.size());
            fitnessEntropyImpactTail = fitnessEntropyImpact.subList(inds.size(),fitnessEntropyImpact.size());
        }

        countingImpactString += state.generation +", " + DoubleListOperations.listToString(countingImpactHead) + ", " + DoubleListOperations.getSum(countingImpactTail) + ", " + DoubleListOperations.listToString(countingImpactTail) + "\n";
        fitnessImpactString += state.generation +", " + DoubleListOperations.listToString(fitnessImpactHead) + ", " + DoubleListOperations.getSum(fitnessImpactTail) + ", " + DoubleListOperations.listToString(fitnessImpactTail) + "\n";
        entropyImpactString += state.generation +", " + DoubleListOperations.listToString(entropyImpactHead) + ", " + DoubleListOperations.getSum(entropyImpactTail) + ", " + DoubleListOperations.listToString(entropyImpactTail) + "\n";
        fitnessEntropyImpactString += state.generation +", " + DoubleListOperations.listToString(fitnessEntropyImpactHead) + ", " + DoubleListOperations.getSum(fitnessEntropyImpactTail) + ", " + DoubleListOperations.listToString(fitnessEntropyImpactTail) + "\n";


        //normalize and print the normalized files
        List<Double> normalizedFitnessImpact = DoubleListOperations.normalizeList(fitnessImpact);
        List<Double> normalizedEntropyImpact = DoubleListOperations.normalizeList(entropyImpact);
        List<Double> normalizedFitnessEntropyImpact = DoubleListOperations.normalizeList(fitnessEntropyImpact);
        List<Double> normalizedFitnessImpactHead = normalizedFitnessImpact.subList(0,inds.size());
        List<Double> normalizedEntropyImpactHead = normalizedEntropyImpact.subList(0,inds.size());
        List<Double> normalizedFitnessEntropyImpactHead = normalizedFitnessEntropyImpact.subList(0,inds.size());
        List<Double> normalizedFitnessImpactTail = new ArrayList<Double>();
        List<Double> normalizedEntropyImpactTail = new ArrayList<Double>();
        List<Double> normalizedFitnessEntropyImpactTail = new ArrayList<Double>();
        if(state.mutationCounter < 0){
            normalizedFitnessImpactTail = normalizedFitnessImpact.subList(normalizedFitnessImpactHead.size(),normalizedFitnessImpact.size());
            normalizedEntropyImpactTail = normalizedEntropyImpact.subList(normalizedEntropyImpactHead.size(),normalizedEntropyImpact.size());
            normalizedFitnessEntropyImpactTail = normalizedFitnessEntropyImpact.subList(normalizedFitnessEntropyImpactHead.size(),normalizedFitnessEntropyImpact.size());
        }
        normalizedFitnessImpactString += state.generation +", " + DoubleListOperations.listToString(normalizedFitnessImpactHead) + ", " + DoubleListOperations.getSum(normalizedFitnessImpactTail) + ", " + DoubleListOperations.listToString(normalizedFitnessImpactTail) + "\n";
        normalizedEntropyImpactString += state.generation +", " + DoubleListOperations.listToString(normalizedEntropyImpactHead) + ", " + DoubleListOperations.getSum(normalizedEntropyImpactTail) + ", " + DoubleListOperations.listToString(normalizedEntropyImpactTail) + "\n";
        normalizedFitnessEntropyImpactString += state.generation +", " + DoubleListOperations.listToString(normalizedFitnessEntropyImpactHead) + ", " + DoubleListOperations.getSum(normalizedFitnessEntropyImpactTail) + ", " + DoubleListOperations.listToString(normalizedFitnessEntropyImpactTail) + "\n";

        //print the entropyFile
        entropyString += state.generation + ", " + IndividualListOperations.calculateEntropy(inds) + "\n";
    }

    @Override
    public void finalStatistics(EvolutionState state, int result) {
        int genomeLength = IndividualListOperations.getGenomeLength(state.population.subpops.get(0).individuals.get(0));
        String headder = "";

        ArrayList<Individual> inds = state.population.subpops.get(0).individuals;

        //print the startingEndingFitness and the startingFitness
        headder = "BestInd,Mean,Median";

        startingEndingFitness = headder +"\n" + startingEndingFitness;
        state.output.println(startingEndingFitness, statLog);

        //fitness and normalizedFitness headder
        headder = "generation";
        for(int i = 1; i <= inds.size(); i++)
            headder += ", " + i;
        fitnessString = headder + "\n" + fitnessString;

        //extend to impact and normalized counterparts headder
        headder += ", m";
        for(int i = -1; i >= state.mutationCounter; i--)
            headder += ", " + i;
        countingImpactString = headder + "\n" + countingImpactString;
        fitnessImpactString = headder + "\n" + fitnessImpactString;
        entropyImpactString = headder + "\n" + entropyImpactString;
        fitnessEntropyImpactString = headder + "\n" + fitnessEntropyImpactString;
        normalizedFitnessImpactString = headder + "\n" + normalizedFitnessImpactString;
        normalizedEntropyImpactString = headder + "\n" + normalizedEntropyImpactString;
        normalizedFitnessEntropyImpactString = headder + "\n" + normalizedFitnessEntropyImpactString;

        //traceID and genome headder
        headder = "generation, indID";
        for(int i = 1; i <= genomeLength; i++)
            headder += ", traceID_" + i;
        traceIDsString = headder + "\n" + traceIDsString;

        //traceID and genome headder
        headder = "generation, indID";
        for(int i = 1; i <= genomeLength; i++)
            headder += ", genome_" + i;
        traceIDsString = headder + "\n" + traceIDsString;
        genomeString = headder + "\n" + genomeString;

        //entropy headder
        headder = "generation, overall entropy, entropy sum";
        for(int i = 1; i <= genomeLength; i++)
            headder += ", entropy ind_" + i;
        entropyString = headder + "\n" + entropyString;


        //print the files
        state.output.println(fitnessString, fitnessLog);

        state.output.println(traceIDsString, traceIDsLog);
        state.output.println(genomeString, genomeLog);

        state.output.println(countingImpactString, countingImpactLog);
        state.output.println(fitnessImpactString, fitnessImpactLog);
        state.output.println(entropyImpactString, entropyImpactLog);
        state.output.println(fitnessEntropyImpactString, fitnessEntropyImpactLog);
        state.output.println(normalizedFitnessImpactString, normalizedFitnessImpactLog);
        state.output.println(normalizedEntropyImpactString, normalizedEntropyImpactLog);
        state.output.println(normalizedFitnessEntropyImpactString, normalizedFitnessEntropyImpactLog);

        state.output.println(entropyString, entropyLog);
    }
}
