package ec.app.TracableVectorProblems.MultiObjective.WFG;

import ec.EvolutionState;
import ec.app.TracableVectorProblems.MultiObjective.MultiObjectiveTraceableDoubleVectorProblem;
import ec.util.Parameter;
import ec.vector.FloatVectorSpecies;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Implements a reference abstract class for all wfg org.uma.test problem
 * Reference: Simon Huband, Luigi Barone, Lyndon While, Phil Hingston
 * A Scalable Multi-objective Test Problem Toolkit.
 * Evolutionary Multi-Criterion Optimization:
 * Third International Conference, EMO 2005.
 * Proceedings, volume 3410 of Lecture Notes in Computer Science
 */

public abstract class WFG extends MultiObjectiveTraceableDoubleVectorProblem {

    public static final String P_K = "K"; //the number of position parameters
    public static final String P_L = "L"; //the number of distance parameters

    /**
     * stores a epsilon default value
     */
    private final float epsilon = (float) 1e-7;

    protected int k; //the number of position parameters
    protected int m; //the number of objective functions
    protected int l; //the number of distance parameters
    protected int[] a;
    protected int[] s;
    protected int d = 1;
    protected Random random = new Random();



    @Override
    public void initialSetup(final EvolutionState state, final Parameter base) {
        //set min and max number of objectives
        this.minNumberObjectives = 2; //the number of objectives for this problem
        this.maxNumberObjectives = 3; //the number of objectives for this problem

        //read K and M, check if they are in the genome size boud
        Parameter p;
        p = new Parameter(P_K);
        if (state.parameters.exists(p, null)) {
            this.k = state.parameters.getInt(p, null, 0);
        } else {
            state.output.fatal("parameter \"K\" not defined. Check your params file.");
        }
        p = new Parameter(P_L);
        if (state.parameters.exists(p, null)) {
            this.l = state.parameters.getInt(p, null, 0);
        } else {
            state.output.fatal("parameter \"L\" not defined. Check your params file.");
        }

        //check if genomeSize and l+k is equal
        int genomeSize = state.parameters.getInt(new Parameter("pop.subpop.0.species.genome-size"), null);
        if(genomeSize != k+l)
            state.output.fatal("The number of distance parameters L="+l+" and the number of position parameters K="+k+" does not match the genomeSize of "+genomeSize);

        //set the bounds dynamically with the population size
        int i = 1;
        while (i < genomeSize) {
            state.parameters.set(new Parameter("pop.subpop.0.species.min-gene."+i), "0");
            state.parameters.set(new Parameter("pop.subpop.0.species.max-gene."+i), 2.0*(i+1)+"");
            i++;
        }
        state.parameters.set(new Parameter("pop.subpop.0.species.min-gene"), "0");
        state.parameters.set(new Parameter("pop.subpop.0.species.max-gene"), 2.0*(i+1)+"");

        //set the number of objectives
        this.m = state.parameters.getInt(new Parameter("multi.fitness.num-objectives"), null);
    }

    /**
     * Gets the x vector
     */
    public float[] calculateX(float[] t) {
        float[] x = new float[m];

        for (int i = 0; i < m - 1; i++) {
            x[i] = Math.max(t[m - 1], a[i]) * (t[i] - (float) 0.5) + (float) 0.5;
        }

        x[m - 1] = t[m - 1];

        return x;
    }

    /**
     * Normalizes a vector (consulte wfg toolkit reference)
     */
    public float[] normalise(double[] z) {
        float[] result = new float[z.length];

        for (int i = 0; i < z.length; i++) {
            float bound = (float) 2.0 * (i + 1);
            result[i] = (float)(z[i] / bound);
            result[i] = (float)correctTo01(result[i]);
        }

        return result;
    }

    /**
     */
    public double correctTo01(double a) {
        double min = (double) 0.0;
        double max = (double) 1.0;

        double minEpsilon = min - epsilon;
        double maxEpsilon = max + epsilon;

        if ((a <= min && a >= minEpsilon) || (a >= min && a <= minEpsilon)) {
            return min;
        } else if ((a >= max && a <= maxEpsilon) || (a <= max && a >= maxEpsilon)) {
            return max;
        } else {
            return a;
        }
    }

    /**
     * Gets a subvector of a given vector
     * (Head inclusive and tail inclusive)
     *
     * @param z the vector
     * @return the subvector
     */
    public float[] subVector(float[] z, int head, int tail) {
        int size = tail - head + 1;
        float[] result = new float[size];

        System.arraycopy(z, head, result, head - head, tail + 1 - head);

        return result;
    }


    //TODO: remove?
    /**
     * Evaluates a solution
     *
     * @param variables The solution to evaluate
     * @return a double [] with the evaluation results
     */
    //abstract public float[] evaluate(float[] variables);
}