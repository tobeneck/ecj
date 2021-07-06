package ec.app.TracableVectorProblems.TracableVectorStatistics;

import ec.EvolutionState;
import ec.Individual;
import ec.Statistics;
import ec.app.TracableVectorProblems.TracableVectorStatistics.ListOperations.IndividualAndGenomeListOperations;
import ec.multiobjective.nsga2.NSGA2MultiObjectiveFitness;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceTuple;
import ec.vector.TracableDataTypes.TraceableString;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static ec.app.TracableVectorProblems.TracableVectorStatistics.ListOperations.IndividualAndGenomeListOperations.*;

public class TraceableVectorStatistics extends Statistics
{
    public static final String P_STARTINGPOP = "starting-population-file";
    public int startingPopulationLog;

    public static final String P_EVALFILE = "eval-file"; //eval file wich contains the information for the python evaluation
    public int evalLog;
    protected String evalString = "generation;individual;fitness;genome\n";



    /**
     * returns the eval string for a given list of genomes
     * @param genomes the genomes to be stringyfied
     * @return string of the input genomes
     */
    protected String getGenomeString(ArrayList<TraceableString> genomes){
        String outString = "";
        for(int i = 0; i < genomes.size(); i++){
            TraceableString currentGene = genomes.get(i);

            //add the value
            outString += currentGene.getValue();

            //add the trace list
            for(int j = 0; j < currentGene.getTraceVector().size(); j++){
                TraceTuple currentTraceTuple = currentGene.getTraceVector().get(j);
                outString += ","+currentTraceTuple.getTraceID()+","+currentTraceTuple.getImpact();
            }

            if( i+1 < genomes.size())
                outString += "|";
        }

        return outString;
    }

    /**
     * returns the eval string for a given list of fitness values
     * @param fitness the fitness to be stringyfied
     * @return string of the input fitness
     */
    protected String getFitnessString(double[] fitness){
        String outString = "";

        for(int j = 0; j < fitness.length; j++){
            outString += fitness[j];
            if(j + 1< fitness.length)
                outString += "|";
        }

        return outString;
    }

    public void setup(final EvolutionState state, final Parameter base) {
        // DO NOT FORGET to call super.setup(...) !!
        super.setup(state,base);

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

        // set up evalFile
        File evalFile = state.parameters.getFile(
                base.push(P_EVALFILE),null);
        if (evalFile!=null) try
        {
            evalLog = state.output.addLog(evalFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    evalFile + ":\n" + i);
        }

    }

    public void postEvaluationStatistics(final EvolutionState state)
    {

        ArrayList<Individual> inds = state.population.subpops.get(0).individuals; //provides fitness and strings to deconstruct
        ArrayList<ArrayList<TraceableString>> genomes = new ArrayList<ArrayList<TraceableString>>(); //provides the genomes (values and traceVectors) of the inds. As the inds have a generic datatype, this is the only way to keep it interchangeable with other traceableDatatypes.
        for(int i = 0; i < inds.size(); i++){
            genomes.add(new ArrayList<TraceableString>());
            String[] genotypeString = getGenotypeString(inds.get(i));
            for(int j = 0; j < genotypeString.length; j++) {
                TraceableString currentInd = new TraceableString();
                currentInd.fromString(genotypeString[j]);
                genomes.get(i).add(currentInd);
            }
        }

        //print the startingPopulationFile and the startingFitness
        if(state.generation == 0){
            // print out the population
            state.population.printPopulation(state,startingPopulationLog);
        }


        //print the eval string
        String generationString = state.generation + "";

        for(int i = 0; i < inds.size(); i++){
            String individualString = i+"";

            String fitnessString = getFitnessString(getFitness(inds.get(i)));

            String genomeString = getGenomeString(genomes.get(i));

            //TODO: add "rank" string, that shows the quality of solution?

            evalString += generationString+";"+individualString+";"+fitnessString+";"+genomeString+"\n";
        }
    }

    @Override
    public void finalStatistics(EvolutionState state, int result) {
        state.output.println(evalString, evalLog);
    }


}
