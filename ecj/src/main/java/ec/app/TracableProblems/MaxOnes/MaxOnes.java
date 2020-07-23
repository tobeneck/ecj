/*
  Copyright 2006 by Sean Luke
  Licensed under the Academic Free License version 3.0
  See the file "LICENSE" for more information
*/


package ec.app.TracableProblems.MaxOnes;
import ec.*;
import ec.simple.*;
import ec.vector.*;

public class MaxOnes extends Problem implements SimpleProblemForm {

    public void evaluate(final EvolutionState state,
                         final Individual ind,
                         final int subpopulation,
                         final int threadnum) {
        if (ind.evaluated) return;

        if (!(ind instanceof TracableBitVectorIndividual))
            state.output.fatal("Whoa!  It's not a (Tracable)BitVectorIndividual!!!", null);


        int sum=0;
        TracableBitVectorIndividual ind2 = (TracableBitVectorIndividual)ind;

        for(int x=0; x<ind2.genome.length; x++)
            sum += (ind2.genome[x].getValue() ? 1 : 0);

        if (!(ind2.fitness instanceof SimpleFitness))
            state.output.fatal("Whoa!  It's not a SimpleFitness!!!",null);
        ((SimpleFitness)ind2.fitness).setFitness(state,
                /// ...the fitness...
                sum/(double)ind2.genome.length,
                ///... is the individual ideal?  Indicate here...
                sum == ind2.genome.length);
        ind2.evaluated = true;

    }
}
