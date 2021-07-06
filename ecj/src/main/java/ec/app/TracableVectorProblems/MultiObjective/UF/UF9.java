package ec.app.TracableVectorProblems.MultiObjective.UF;

import ec.EvolutionState;
import ec.app.TracableVectorProblems.MultiObjective.MultiObjectiveTraceableDoubleVectorProblem;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;

public class UF9 extends MultiObjectiveTraceableDoubleVectorProblem {

    public static final String P_EPSILON = "epsilon";
    double epsilon;

    public void initialSetup(final EvolutionState state, final Parameter base) {
        this.minNumberObjectives = 3; //the number of objectives for this problem
        this.maxNumberObjectives = 3; //the number of objectives for this problem
        this.problemName = "UF9"; //the name of the problem

        Parameter p;
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

        int count1, count2, count3;
        double sum1, sum2, sum3, yj;
        sum1   = sum2 = sum3 = 0.0;
        count1 = count2 = count3 = 0;

        for (int j = 3 ; j <= genomeLength; j++) {
            yj = x[j-1] - 2.0*x[1]*Math.sin(2.0*Math.PI*x[0]+j*Math.PI/genomeLength);
            yj = x[j-1] - 2.0*x[1]*Math.sin(2.0*Math.PI*x[0]+j*Math.PI/genomeLength);
            if(j % 3 == 1) {
                sum1  += yj*yj;
                count1++;
            } else if(j % 3 == 2) {
                sum2  += yj*yj;
                count2++;
            } else {
                sum3  += yj*yj;
                count3++;
            }
        }

        yj = (1.0+epsilon)*(1.0-4.0*(2.0*x[0]-1.0)*(2.0*x[0]-1.0));
        if (yj < 0.0)
            yj = 0.0;

        f[0] = 0.5*(yj + 2*x[0])*x[1]		+ 2.0*sum1 / (double)count1;
        f[1] = 0.5*(yj - 2*x[0] + 2.0)*x[1] + 2.0*sum2 / (double)count2;
        f[2] = 1.0 - x[1]                   + 2.0*sum3 / (double)count3 ;


        return f;
    }
}
