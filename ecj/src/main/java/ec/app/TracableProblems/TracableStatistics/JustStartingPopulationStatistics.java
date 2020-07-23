package ec.app.TracableProblems.TracableStatistics;

import ec.EvolutionState;
import ec.Individual;
import ec.Statistics;
import ec.app.TracableProblems.TracableStatistics.ListOperations.DoubleListOperations;
import ec.app.TracableProblems.TracableStatistics.ListOperations.IndividualListOperations;
import ec.util.Parameter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class JustStartingPopulationStatistics extends Statistics {
    // The parameter string and log number of the file for our readable population
    public static final String P_STATFILE = "stat-file"; //TODO: rename
    public int statLog;

    public static final String P_STARTINGPOP = "starting-population-file";
    public int startingPopulationLog;

    String startingEndingFitness = "";


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
    }



    public void postEvaluationStatistics(final EvolutionState state)
    {
        // be certain to call the hook on super!
        super.postEvaluationStatistics(state);

        ArrayList<Individual> inds = state.population.subpops.get(0).individuals;

        //print the startingPopulationFile and the startingFitness
        if(state.generation == 0){
            // print out the population
            state.population.subpops.get(0).printSubpopulation(state,startingPopulationLog);
            startingEndingFitness = "" + IndividualListOperations.getBest(inds) + "," + IndividualListOperations.getMean(inds) +","+IndividualListOperations.getMedian(inds);
        }
    }

    @Override
    public void finalStatistics(EvolutionState state, int result) {
        int genomeLength = IndividualListOperations.getGenomeLength(state.population.subpops.get(0).individuals.get(0));
        String headder = "BestInd,Mean,Median";

        startingEndingFitness = headder +"\n" + startingEndingFitness;

        state.output.println(startingEndingFitness, statLog);
    }
}
