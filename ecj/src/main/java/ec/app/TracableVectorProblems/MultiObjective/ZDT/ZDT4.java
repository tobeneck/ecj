package ec.app.TracableVectorProblems.MultiObjective.ZDT;

import ec.EvolutionState;
import ec.app.TracableVectorProblems.MultiObjective.ZDT.ZDT1;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;

public class ZDT4 extends ZDT1
{
    public void setup(final EvolutionState state, final Parameter base) { }


    /**
     * Returns the value of the ZDT4 function G.
     *
     * @param genome Genome
     */
    @Override
    public double evalG(TraceableDouble[] genome) {
        double g = 0.0;
        for (int var = 1; var < genome.length; var++) {
            g += Math.pow(genome[var].getValue(), 2.0) +
                    -10.0 * Math.cos(4.0 * Math.PI * genome[var].getValue());
        }

        double constant = 1.0 + 10.0 * (genome.length - 1); //TODO: should this be 91?
        return g + constant;
    }

    /**
     * Returns the value of the ZDT4 function H.
     *
     * @param f1 First argument of the function H.
     * @param g Second argument of the function H.
     */
    @Override
    public double evalH(double f1, double g) {
        return 1.0 - Math.sqrt(f1 / g);
    }

    //evaluate function stays the same
}