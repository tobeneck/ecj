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

public class RosenbrockFunction extends Problem implements SimpleProblemForm
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
        double value = 0;

        // Compute the Rosenbrock function for our genome
//        for( int i = 1 ; i < len ; i++ )
//            value += 100*(genome[i-1].getValue()*genome[i-1].getValue()-genome[i].getValue())*
//                (genome[i-1].getValue()*genome[i-1].getValue()-genome[i].getValue()) +
//                (1-genome[i-1].getValue())*(1-genome[i-1].getValue());

        //Wikipedia variante: a=1, b=100
        for( int i = 0 ; i < len - 1 ; i++ )
            value += 100*(genome[i+1].getValue()-genome[i].getValue()*genome[i].getValue())*
                    (genome[i+1].getValue()-genome[i].getValue()*genome[i].getValue()) +
                    (1-genome[i].getValue())*(1-genome[i].getValue());


        // Rosenbrock is a minimizing function which does not drop below 0. 
        // But SimpleFitness requires a maximizing function -- where 0 is worst
        // and 1 is best.  To use SimpleFitness, we must convert the function.
        // This is the Koza style of doing it:

        //global minimum f(x) = 0. Fitness = 0 - 1, higher is better
        value = 1.0 / ( 1.0 + value );
        ((SimpleFitness)(ind.fitness)).setFitness( state, value, value==1.0 );
    
        ind.evaluated = true;
        }
    }
