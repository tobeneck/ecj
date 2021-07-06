package ec.app.TracableVectorProblems.MultiObjective.UF;

import ec.EvolutionState;
import ec.app.TracableVectorProblems.MultiObjective.MultiObjectiveTraceableDoubleVectorProblem;
import ec.multiobjective.MultiObjectiveFitness;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;

public class UF1 extends MultiObjectiveTraceableDoubleVectorProblem {

    public void initialSetup(final EvolutionState state, final Parameter base) {
        this.minNumberObjectives = 2; //the number of objectives for this problem
        this.maxNumberObjectives = 2; //the number of objectives for this problem
        this.problemName = "UF1"; //the name of the problem
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

        int genomeLength = genome.length;

        double[] f = new double[numberObjectives]; //the multi objective fitness

        double[] x = new double[genome.length];
        for (int i = 0; i < genomeLength; i++) {
            x[i] = genome[i].getValue();
        }

        int count1, count2;
        double sum1, sum2, yj;
        sum1 = sum2 = 0.0;
        count1 = count2 = 0;

        for (int j = 2; j <= genomeLength; j++) {
            yj = x[j - 1] - Math.sin(6.0 * Math.PI * x[0] + j * Math.PI / genomeLength);
            yj = yj * yj;
            if (j % 2 == 0) {
                sum2 += yj;
                count2++;
            } else {
                sum1 += yj;
                count1++;
            }
        }

        f[0] = x[0] + 2.0 * sum1 / (double) count1;
        f[1] = 1.0 - Math.sqrt(x[0]) + 2.0 * sum2 / (double) count2;

        return f;

    }
}
