package ec.app.TracableVectorProblems.MultiObjective.ZDT;

import ec.EvolutionState;
import ec.util.Parameter;

public class ZDT2 extends ZDT1
{
    public void initialSetup(final EvolutionState state, final Parameter base) {
        super.initialSetup(state, base);
        this.minNumberObjectives = 2; //the number of objectives for this problem
        this.maxNumberObjectives = 2; //the number of objectives for this problem
        this.problemName = "ZDT2"; //the name of the problem
    }

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
