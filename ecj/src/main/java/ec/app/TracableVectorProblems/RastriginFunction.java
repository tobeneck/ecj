/*
  Copyright 2006 by Sean Luke
  Licensed under the Academic Free License version 3.0
  See the file "LICENSE" for more information
*/


package ec.app.TracableVectorProblems;

import ec.EvolutionState;
import ec.Individual;
import ec.Problem;
import ec.simple.SimpleFitness;
import ec.simple.SimpleProblemForm;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;
import ec.vector.TraceableDoubleVectorIndividual;

public class RastriginFunction extends Problem implements SimpleProblemForm
{
    public void setup(final EvolutionState state, final Parameter base) { }

    public void evaluate(final EvolutionState state,
                         final Individual ind,
                         final int subpopulation,
                         final int threadnum)
    {
        if( !( ind instanceof TraceableDoubleVectorIndividual ) )
            state.output.fatal( "The individuals for this problem should be TraceableDoubleVectorIndividuals." );

        TraceableDouble[] genome = ((TraceableDoubleVectorIndividual)ind).genome;
        int len = genome.length;
        int A = 10;
        double value = (double) (A * len);


        double sum = 0;
        //Wikipedia variante:
        for( int i = 0 ; i < len - 1 ; i++ )
            sum += 100*(genome[i].getValue()*genome[i].getValue()) - A*Math.cos(2*Math.PI*genome[i].getValue());

        value = value + sum;

        //global minimum f(x) = 0. Fitness = 0 - 1, higher is better
        value = 1.0 / ( 1.0 + value );
        ((SimpleFitness)(ind.fitness)).setFitness( state, value, value==1.0 ); //TODO: is this right?

        ind.evaluated = true;
    }
}
