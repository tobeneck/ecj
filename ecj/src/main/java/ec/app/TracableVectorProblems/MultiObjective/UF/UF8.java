package ec.app.TracableVectorProblems.MultiObjective.UF;

import ec.EvolutionState;
import ec.app.TracableVectorProblems.MultiObjective.MultiObjectiveTraceableDoubleVectorProblem;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;

public class UF8 extends MultiObjectiveTraceableDoubleVectorProblem {

    public void initialSetup(final EvolutionState state, final Parameter base) {
        this.minNumberObjectives = 3; //the number of objectives for this problem
        this.maxNumberObjectives = 3; //the number of objectives for this problem
        this.problemName = "UF8"; //the name of the problem
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
        sum1 = sum2 = sum3 = 0.0;
        count1 = count2 = count3 = 0;

        for (int j = 3; j <= genomeLength; j++) {
            yj =
                    x[j - 1]
                            - 2.0 * x[1] * Math.sin(2.0 * Math.PI * x[0] + j * Math.PI / genomeLength);
            if (j % 3 == 1) {
                sum1 += yj * yj;
                count1++;
            } else if (j % 3 == 2) {
                sum2 += yj * yj;
                count2++;
            } else {
                sum3 += yj * yj;
                count3++;
            }
        }

        f[0] = Math.cos(0.5 * Math.PI * x[0]) * Math.cos(0.5 * Math.PI * x[1]) + 2.0 * sum1 / (double) count1;
        f[1] = Math.cos(0.5 * Math.PI * x[0]) * Math.sin(0.5 * Math.PI * x[1]) + 2.0 * sum2 / (double) count2;
        f[2] = Math.sin(0.5 * Math.PI * x[0]) + 2.0 * sum3 / (double) count3;


        return f;
    }
}
