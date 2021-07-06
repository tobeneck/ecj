package ec.app.TracableVectorProblems.MultiObjective.ZDT;

import ec.EvolutionState;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;

public class ZDT6 extends ZDT1
{
    public void initialSetup(final EvolutionState state, final Parameter base) {
        super.initialSetup(state, base);
        this.minNumberObjectives = 2; //the number of objectives for this problem
        this.maxNumberObjectives = 2; //the number of objectives for this problem
        this.problemName = "ZDT6"; //the name of the problem
    }

    /**
     * Returns the value of the ZDT6 function G.
     *
     * @param genome Genome
     */
    protected double evalG(TraceableDouble[] genome) {
        double g = 0.0;
        for (int var = 1; var < genome.length; var++) {
            g += genome[var].getValue();
        }
        g = g / (genome.length - 1);
        g = Math.pow(g, 0.25);
        g = 9.0 * g;
        g = 1.0 + g;
        return g;
    }


    /**
     * Returns the value of the ZDT6 function H.
     *
     * @param f1 First argument of the function H.
     * @param g Second argument of the function H.
     */
    protected double evalH(double f1, double g) {
        return 1.0 - Math.pow((f1 / g), 2.0);
    }

    /**
     * returns the fitness of the current problem. Mainly combines the output of evalF and evalH
     * @param genome the genome to be evaluated
     * @param numberObjectives the number of Objectives in the problem
     * @param state the evolutionary state for error output
     * @return the fitness (double array of size 2) of the genome
     */
    @Override
    protected double[] getFitness(TraceableDouble[] genome, int numberObjectives, EvolutionState state){

        double[] f = new double[numberObjectives]; //the multi objective fitness

        f[0] = 1 - Math.exp(-4 * genome[0].getValue()) * Math.pow(Math.sin(6 * Math.PI * genome[0].getValue()), 6);
        double g = this.evalG(genome);
        double h = this.evalH(f[0], g);
        f[1] = h * g;

        return f;

    }
}
