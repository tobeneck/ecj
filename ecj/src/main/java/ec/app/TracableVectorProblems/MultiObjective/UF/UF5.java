package ec.app.TracableVectorProblems.MultiObjective.UF;

import ec.EvolutionState;
import ec.app.TracableVectorProblems.MultiObjective.MultiObjectiveTraceableDoubleVectorProblem;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;

public class UF5 extends MultiObjectiveTraceableDoubleVectorProblem {

    public static final String P_N = "N";
    int n;

    public static final String P_EPSILON = "epsilon";
    double epsilon;

    public void initialSetup(final EvolutionState state, final Parameter base) {
        this.minNumberObjectives = 2; //the number of objectives for this problem
        this.maxNumberObjectives = 2; //the number of objectives for this problem
        this.problemName = "UF5"; //the name of the problem

        //read the other external parameters
        Parameter p;

        p = new Parameter(P_N);
        if (state.parameters.exists(p, null)) {
            this.n = state.parameters.getInt(p, null, 0);
        } else {
            state.output.fatal("parameter \"N\" not defined. Check your params file.");
        }
        p = new Parameter(P_EPSILON);
        if (state.parameters.exists(p, null)) {
            this.epsilon = state.parameters.getDouble(p, null, 0);
        } else {
            state.output.fatal("parameter \"epsilon\" not defined. Check your params file.");
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
    protected double[] getFitness(TraceableDouble[] genome, int numberObjectives, EvolutionState state){

        double[] f = new double[numberObjectives]; //the multi objective fitness

        int genomeLength = genome.length;

        double[] x = new double[genome.length];
        for (int i = 0; i < genomeLength; i++) {
            x[i] = genome[i].getValue();
        }

        int count1, count2;
        double sum1, sum2, yj, hj;
        sum1 = sum2 = 0.0;
        count1 = count2 = 0;

        for (int j = 2; j <= genomeLength; j++) {
            yj = x[j - 1] - Math.sin(6.0 * Math.PI * x[0] + j * Math.PI / genomeLength);
            hj = Math.abs(yj) / (1.0 + Math.exp(2.0 * Math.abs(yj)));
            if (j % 2 == 0) {
                sum2 += hj;
                count2++;
            } else {
                sum1 += hj;
                count1++;
            }
        }

        f[0] = x[0] + 2.0 * sum1 / (double) count1;
        f[1] = 1.0 - x[0] * x[0] + 2.0 * sum2 / (double) count2;

        return f;
    }
}
