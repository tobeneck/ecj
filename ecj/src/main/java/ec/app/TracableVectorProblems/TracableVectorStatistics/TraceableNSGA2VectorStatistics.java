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

public class TraceableNSGA2VectorStatistics extends TraceableVectorStatistics
{
    @Override
    public void setup(final EvolutionState state, final Parameter base)
    {
        //call super to initialize the logs and files correctly
        super.setup(state, base);

        //update hte headline of the eval string
        evalString = "generation;individual;fitness;genome;rank;sparsity\n";
    }

    @Override
    public void postEvaluationStatistics(final EvolutionState state) //this is basically the same method as the parent, but with rank and sparsity automatically atattched to the file
    {
        // be certain to call the hook on super!
        //super.postEvaluationStatistics(state); //actually, don't call super because we want the eval string to be updated there!

        ArrayList<Individual> inds = state.population.subpops.get(0).individuals; //provides fitness and strings to deconstruct

        if(!(inds.get(0).fitness instanceof NSGA2MultiObjectiveFitness))
            state.output.fatal("Fitness of the Individual needs to be an instance of NSGA2MultiObjectiveFitness, please check your config files!\n");


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

        //print the startingPopulationFile
        if(state.generation == 0){
            state.population.printPopulation(state,startingPopulationLog);
        }



        //print the eval string
        String generationString = state.generation + "";

        for(int i = 0; i < inds.size(); i++){
            String individualString = i+"";

            String fitnessString = getFitnessString(getFitness(inds.get(i)));

            String genomeString = getGenomeString(genomes.get(i));

            String rankString = ((NSGA2MultiObjectiveFitness)inds.get(i).fitness).rank+"";

            String sparsityString = ((NSGA2MultiObjectiveFitness)inds.get(i).fitness).sparsity+"";

            evalString += generationString+";"+individualString+";"+fitnessString+";"+genomeString+";"+rankString+";"+sparsityString+"\n";
        }
    }


}
