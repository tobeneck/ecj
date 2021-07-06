package ec.app.TracableVectorProblems;

        import ec.EvolutionState;
        import ec.Individual;
        import ec.Problem;
        import ec.multiobjective.nsga2.NSGA2MultiObjectiveFitness;
        import ec.simple.SimpleFitness;
        import ec.simple.SimpleProblemForm;
        import ec.util.Parameter;
        import ec.vector.TracableDataTypes.TraceableDouble;
        import ec.vector.TraceableDoubleVectorIndividual;

public class SchafferFunctionN1 extends Problem implements SimpleProblemForm
{
    public void setup(final EvolutionState state, final Parameter base) { }

    public void evaluate(final EvolutionState state,
                         final Individual ind,
                         final int subpopulation,
                         final int threadnum)
    {
        if( !( ind instanceof TraceableDoubleVectorIndividual) )
            state.output.fatal( "The individuals for this problem should be TraceableDoubleVectorIndividuals." );


        TraceableDouble[] genome = ((TraceableDoubleVectorIndividual)ind).genome;
        int len = genome.length;

        if( len != 1 )
            state.output.fatal( "The SchafferFunctionN1 only works with a genome length of 1, not of: "+ len );

        double x = genome[0].getValue();
        double obj1 = x*x;
        double obj2 = (x-2)*(x-2);


        // Rosenbrock is a minimizing function which does not drop below 0.
        // But SimpleFitness requires a maximizing function -- where 0 is worst
        // and 1 is best.  To use SimpleFitness, we must convert the function.
        // This is the Koza style of doing it:

        //global minimum f(x) = 0. Fitness = 0 - 1, higher is better
        obj1 = 1.0 / ( 1.0 + obj1 );
        obj2 = 1.0 / ( 1.0 + obj2 );


        double[] newObjective = {obj1, obj2};
        ((NSGA2MultiObjectiveFitness) (ind.fitness)).setObjectives(state, newObjective);

        ind.evaluated = true;
    }
}
