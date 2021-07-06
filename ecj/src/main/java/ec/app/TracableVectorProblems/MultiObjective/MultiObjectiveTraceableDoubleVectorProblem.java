package ec.app.TracableVectorProblems.MultiObjective;

import ec.EvolutionState;
import ec.Individual;
import ec.Problem;
import ec.multiobjective.MultiObjectiveFitness;
import ec.simple.SimpleProblemForm;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;
import ec.vector.TraceableDoubleVectorIndividual;

public abstract class MultiObjectiveTraceableDoubleVectorProblem extends Problem implements SimpleProblemForm
{
    protected int minGenomeLength = 1; //the min genome length for the problem to check against
    protected int maxGenomeLength = Integer.MAX_VALUE; //the max genome length for the problem to check against

    protected int minNumberObjectives = 2; //the min number of objectives for this problem
    protected int maxNumberObjectives = 2; //the max number of objectives for this problem

    protected String problemName = "dummy problem"; //the name of the problem

    /**
     * initial setup only called in the first generation
     *
     */
    abstract protected void initialSetup(final EvolutionState state, final Parameter base);

    public void setup(final EvolutionState state, final Parameter base) {
        super.setup(state, base);

        if(state.generation == 0)
            initialSetup(state, base);

    }


    /**
     * returns the fitness of the current problem. Mainly combines the output of evalF and evalH
     * @param genome the genome to be evaluated
     * @param numberObjectives the number of Objectives in the problem
     * @param state the evolutionary state for error output
     * @return the fitness (double array of size 2) of the genome
     */
    protected double[] getFitness(TraceableDouble[] genome, int numberObjectives, EvolutionState state){
        double f[] = new double[0];

        //please override this in the child file

        return f;
    }

    public void evaluate(final EvolutionState state,
                         final Individual ind,
                         final int subpopulation,
                         final int threadnum)
    {
        if( !( ind instanceof TraceableDoubleVectorIndividual ) )
            state.output.fatal( "The individuals for this problem should be TraceableDoubleVectorIndividuals." );



        TraceableDouble[] genome = ((TraceableDoubleVectorIndividual)ind).genome;

        //check the bounds of the problem
        int genomeLength = genome.length;
        int numberObjectives = ((MultiObjectiveFitness) (ind.fitness)).getNumObjectives();
        if( genomeLength < minGenomeLength || genomeLength > maxGenomeLength )
            state.output.fatal( "The "+problemName+" problem only works with a genome length of "+minGenomeLength+" to "+maxGenomeLength+", not of: "+ genome.length );
        if( numberObjectives < minNumberObjectives || numberObjectives > maxNumberObjectives)
            state.output.fatal( "The "+problemName+" problem only works with "+minNumberObjectives+" to "+maxNumberObjectives+" objectives, not with: "+ ((MultiObjectiveFitness) (ind.fitness)).getNumObjectives() );

        //set the new fitness values from the getFitness function
        ((MultiObjectiveFitness) (ind.fitness)).setObjectives(state, getFitness(genome, numberObjectives, state));

        double[] results = getFitness(genome, numberObjectives, state);

        ind.evaluated = true;
    }
}
