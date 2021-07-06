package ec.app.TracableVectorProblems.MultiObjective.WFG;

import ec.EvolutionState;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;

public class WFG1 extends WFG {

    @Override
    public void initialSetup(final EvolutionState state, final Parameter base) {
        super.initialSetup(state, base);
        this.problemName = "WFG1"; //the name of the problem

        //set up s and a
        this.s = new int[m];
        for (int i = 0; i < m; i++) {
            s[i] = 2 * (i + 1);
        }

        this.a = new int[m - 1];
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
        y = t3(y, state);

        y = t4(y, k, m);


        double[] result = new double[m];
        float[] x = calculateX(y);
        for (int m = 1; m <= this.m - 1; m++) {
            result[m - 1] = d * x[this.m - 1] + s[m - 1] * (new Shapes()).convex(x, m);
        }

        result[m - 1] = d * x[m - 1] + s[m - 1] * (new Shapes()).mixed(x, 5, (float) 1.0);



        return result;

    }

    /**
     * WFG1 t1 transformation
     */
    public float[] t1(float[] z, int k) {
        float[] result = new float[z.length];

        System.arraycopy(z, 0, result, 0, k);

        for (int i = k; i < z.length; i++) {
            result[i] = (new Transformations()).sLinear(z[i], (float) 0.35);
        }

        return result;
    }

    /**
     * WFG1 t2 transformation
     */
    public float[] t2(float[] z, int k) {
        float[] result = new float[z.length];

        System.arraycopy(z, 0, result, 0, k);

        for (int i = k; i < z.length; i++) {
            result[i] = (new Transformations()).bFlat(z[i], (float) 0.8, (float) 0.75, (float) 0.85);
        }

        return result;
    }

    /**
     * WFG1 t3 transformation
     *
     */
    public float[] t3(float[] z, EvolutionState state){
        float[] result = new float[z.length];

        for (int i = 0; i < z.length; i++) {
            result[i] = (new Transformations()).bPoly(z[i], (float) 0.02, state);
        }

        return result;
    }

    /**
     * WFG1 t4 transformation
     */
    public float[] t4(float[] z, int k, int M) {
        float[] result = new float[M];
        float[] w = new float[z.length];

        for (int i = 0; i < z.length; i++) {
            w[i] = (float) 2.0 * (i + 1);
        }

        for (int i = 1; i <= M - 1; i++) {
            int head = (i - 1) * k / (M - 1) + 1;
            int tail = i * k / (M - 1);
            float[] subZ = subVector(z, head - 1, tail - 1);
            float[] subW = subVector(w, head - 1, tail - 1);

            result[i - 1] = (new Transformations()).rSum(subZ, subW);
        }

        int head = k + 1 - 1;
        int tail = z.length - 1;
        float[] subZ = subVector(z, head, tail);
        float[] subW = subVector(w, head, tail);
        result[M - 1] = (new Transformations()).rSum(subZ, subW);

        return result;
    }

}
