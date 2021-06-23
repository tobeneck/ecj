package ec.app.TracableVectorProblems.MultiObjective.ZDT;

import ec.EvolutionState;
import ec.Individual;
import ec.Problem;
import ec.multiobjective.nsga2.NSGA2MultiObjectiveFitness;
import ec.simple.SimpleProblemForm;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;
import ec.vector.TraceableDoubleVectorIndividual;

public class ZDT6 extends Problem implements SimpleProblemForm
{
    public void setup(final EvolutionState state, final Parameter base) { }

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

    public void evaluate(final EvolutionState state,
                         final Individual ind,
                         final int subpopulation,
                         final int threadnum)
    {
        if( !( ind instanceof TraceableDoubleVectorIndividual ) )
            state.output.fatal( "The individuals for this problem should be TraceableDoubleVectorIndividuals." );


        TraceableDouble[] genome = ((TraceableDoubleVectorIndividual)ind).genome;
        int len = genome.length;

        if( len != 2 )
            state.output.fatal( "The ZDT6 problem only works with a genome length of 2, not of: "+ len );



        double[] f = new double[len];

        f[0] = 1 - Math.exp(-4 * genome[0].getValue()) * Math.pow(Math.sin(6 * Math.PI * genome[0].getValue()), 6);
        double g = this.evalG(genome);
        double h = this.evalH(f[0], g);
        f[1] = h * g;

        double obj1 = f[0];
        double obj2 = f[1];

        double[] newObjective = {obj1, obj2};
        ((NSGA2MultiObjectiveFitness) (ind.fitness)).setObjectives(state, newObjective);

        ind.evaluated = true;
    }
}
