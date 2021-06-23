package ec.app.TracableVectorProblems.MultiObjective.ZDT;

import ec.EvolutionState;
import ec.Individual;
import ec.Problem;
import ec.multiobjective.nsga2.NSGA2MultiObjectiveFitness;
import ec.simple.SimpleProblemForm;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;
import ec.vector.TraceableDoubleVectorIndividual;

public class ZDT3 extends ZDT1
{
    public void setup(final EvolutionState state, final Parameter base) { }


    /**
     * Returns the value of the ZDT3 function H.
     *
     * @param f1 First argument of the function H.
     * @param g Second argument of the function H.
     */
    @Override
    protected double evalH(double f1, double g) {
        double h ;
        h = 1.0 - Math.sqrt(f1 / g)
                - (f1 / g) * Math.sin(10.0 * Math.PI * f1);
        return h;
    }

    //evaluate function stays the same
}
