package ec.app.TracableVectorProblems.MultiObjective.ZDT;

import ec.EvolutionState;
import ec.util.Parameter;

public class ZDT3 extends ZDT1
{
    public void initialSetup(final EvolutionState state, final Parameter base) {
        super.initialSetup(state, base);
        this.minNumberObjectives = 2; //the number of objectives for this problem
        this.maxNumberObjectives = 2; //the number of objectives for this problem
        this.problemName = "ZDT3"; //the name of the problem
    }

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
