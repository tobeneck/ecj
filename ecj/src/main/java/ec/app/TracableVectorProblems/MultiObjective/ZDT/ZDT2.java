package ec.app.TracableVectorProblems.MultiObjective.ZDT;

import ec.EvolutionState;
import ec.Individual;
import ec.Problem;
import ec.multiobjective.nsga2.NSGA2MultiObjectiveFitness;
import ec.simple.SimpleProblemForm;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;
import ec.vector.TraceableDoubleVectorIndividual;

public class ZDT2 extends ZDT1
{
    public void setup(final EvolutionState state, final Parameter base) { }


    /**
     * Returns the value of the ZDT2 function H.
     *
     * @param f1 First argument of the function H.
     * @param g Second argument of the function H.
     */
    @Override
    public double evalH(double f1, double g) {
        return 1.0 - Math.pow(f1 / g, 2.0);
    }


    //evaluate function stays the same
}
