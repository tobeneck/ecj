package ec.app.TracableVectorProblems.MultiObjective.WFG;

import ec.EvolutionState;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;

public class WFG9 extends WFG {

    @Override
    public void initialSetup(final EvolutionState state, final Parameter base) {
        super.initialSetup(state, base);
        this.problemName = "WFG9"; //the name of the problem

        //set up s and a
        s = new int[m];
        for (int i = 0; i < m; i++) {
            s[i] = 2 * (i + 1);
        }

        a = new int[m - 1];
        for (int i = 0; i < m - 1; i++) {
            a[i] = 1;
        }
    }


    /**
     * returns the fitness of the current problem. Mainly combines the output of evalF and evalH
     * @param genome the genome to be evaluated
     * @param numberObjectives the number of Objectives in the problem
     * @param state the evolutionary state for error output
     * @return the fitness (double array of size 2) of the genome
     */
    @Override
    protected double[] getFitness(TraceableDouble[] genome, int numberObjectives, EvolutionState state){ //numberObjectives is not really needed here

        int genomeLength = genome.length;

        double[] z = new double[genomeLength]; //the genome, but only the double values without the traceID

        int k = genomeLength - numberObjectives + 1;

        for (int i = 0; i < genomeLength; i++) {
            z[i] = genome[i].getValue();
        }



        float[] y = normalise(z);
        y = t1(y, k);
        y = t2(y, k);
        y = t3(y, k, m);

        double[] result = new double[m];
        float[] x = calculateX(y);
        for (int m = 1; m <= this.m; m++) {
            result[m - 1] = d * x[this.m - 1] + s[m - 1] * (new Shapes()).concave(x, m);
        }

        return result;
    }

    /**
     * WFG9 t1 transformation
     */
    public float[] t1(float[] z, int k) {
        float[] result = new float[z.length];
        float[] w = new float[z.length];

        for (int i = 0; i < w.length; i++) {
            w[i] = (float) 1.0;
        }

        for (int i = 0; i < z.length - 1; i++) {
            int head = i + 1;
            int tail = z.length - 1;
            float[] subZ = subVector(z, head, tail);
            float[] subW = subVector(w, head, tail);
            float aux = (new Transformations()).rSum(subZ, subW);
            result[i] = (new Transformations())
                    .bParam(z[i], aux, (float) 0.98 / (float) 49.98, (float) 0.02, (float) 50);
        }

        result[z.length - 1] = z[z.length - 1];
        return result;
    }

    /**
     * WFG9 t2 transformation
     */
    public float[] t2(float[] z, int k) {
        float[] result = new float[z.length];

        for (int i = 0; i < k; i++) {
            result[i] = (new Transformations()).sDecept(z[i], (float) 0.35, (float) 0.001, (float) 0.05);
        }

        for (int i = k; i < z.length; i++) {
            result[i] = (new Transformations()).sMulti(z[i], 30, 95, (float) 0.35);
        }

        return result;
    }

    /**
     * WFG9 t3 transformation
     */
    public float[] t3(float[] z, int k, int M) {
        float[] result = new float[M];

        for (int i = 1; i <= M - 1; i++) {
            int head = (i - 1) * k / (M - 1) + 1;
            int tail = i * k / (M - 1);
            float[] subZ = subVector(z, head - 1, tail - 1);
            result[i - 1] = (new Transformations()).rNonsep(subZ, k / (M - 1));
        }

        int head = k + 1;
        int tail = z.length;
        int l = z.length - k;
        float[] subZ = subVector(z, head - 1, tail - 1);
        result[M - 1] = (new Transformations()).rNonsep(subZ, l);

        return result;
    }
}
