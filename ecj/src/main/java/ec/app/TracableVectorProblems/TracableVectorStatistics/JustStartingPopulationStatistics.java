package ec.app.TracableVectorProblems.TracableVectorStatistics;

import ec.EvolutionState;
import ec.Individual;
import ec.Statistics;
import ec.app.TracableVectorProblems.TracableVectorStatistics.ListOperations.IndividualAndGenomeListOperations;
import ec.util.Parameter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class JustStartingPopulationStatistics extends Statistics {
    // The parameter string and log number of the file for our readable population
    public static final String P_STATFILE = "stat-file"; //TODO: rename
    public int statLog;

    public static final String P_STARTINGPOP = "starting-population-file";
    public int startingPopulationLog;

    public static final String P_STARTINGFITNESSFILE = "starting-fitness-file";
    public int startingFitnessLog;

    String startingEndingFitness = "";
    String startingFitness = "";


    public void setup(final EvolutionState state, final Parameter base)
    {
        // DO NOT FORGET to call super.setup(...) !!
        super.setup(state,base);

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
                base.push(P_STARTINGFITNESSFILE),null);
        if (fitnessFile!=null) try
        {
            startingFitnessLog = state.output.addLog(fitnessFile,true);
        }
        catch (IOException i)
        {
            state.output.fatal("An IOException occurred while trying to create the log " +
                    fitnessFile + ":\n" + i);
        }
    }



    public void postEvaluationStatistics(final EvolutionState state)
    {
        // be certain to call the hook on super!
        super.postEvaluationStatistics(state);

        ArrayList<Individual> inds = state.population.subpops.get(0).individuals;

        //print the startingPopulationFile, the startingEndingFitness and the startingFitness
        if(state.generation == 0){
            // print out the population
            state.population.subpops.get(0).printSubpopulation(state,startingPopulationLog);
            startingEndingFitness = "" + IndividualAndGenomeListOperations.getBest(inds) + "," + IndividualAndGenomeListOperations.getMean(inds) +","+ IndividualAndGenomeListOperations.getMedian(inds);

            //print the fitnessFile
            String fitnessOut = state.generation +"";
            for(int i = 0; i < inds.size(); i++) {
                fitnessOut += "," + inds.get(i).fitness.fitness();
            }
            startingFitness += fitnessOut;
        }
    }

    @Override
    public void finalStatistics(EvolutionState state, int result) {
        int genomeLength = IndividualAndGenomeListOperations.getGenomeLength(state.population.subpops.get(0).individuals.get(0));
        String headder = "BestInd,Mean,Median";

        startingEndingFitness = headder +"\n" + startingEndingFitness;

        state.output.println(startingEndingFitness, statLog);

        //fitness headder
        ArrayList<Individual> inds = state.population.subpops.get(0).individuals;
        headder = "generation";
        for(int i = 1; i <= inds.size(); i++)
            headder += ", " + i;
        startingFitness = headder + "\n" + startingFitness;

        state.output.println(startingFitness, startingFitnessLog);
    }
}
