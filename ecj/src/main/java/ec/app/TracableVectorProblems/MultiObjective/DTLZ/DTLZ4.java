package ec.app.TracableVectorProblems.MultiObjective.DTLZ;

import ec.EvolutionState;
import ec.app.TracableVectorProblems.MultiObjective.MultiObjectiveTraceableDoubleVectorProblem;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;

public class DTLZ4 extends MultiObjectiveTraceableDoubleVectorProblem {

    public void initialSetup(final EvolutionState state, final Parameter base) {
        this.minNumberObjectives = 2; //the number of objectives for this problem
        this.maxNumberObjectives = 3; //the number of objectives for this problem
        this.problemName = "DTLZ4"; //the name of the problem
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

        double[] x = new double[genomeLength] ;

        for (int i = 0; i < genomeLength; i++) {
            x[i] = genome[i].getValue();
        }

        int k = genomeLength - numberObjectives + 1;



        double alpha = 100.0;

        double g = 0.0;
        for (int i = genomeLength - k; i < genomeLength; i++) {
            g += (x[i] - 0.5) * (x[i] - 0.5);
        }

        for (int i = 0; i < numberObjectives; i++) {
            f[i] = 1.0 + g;
        }

        for (int i = 0; i < numberObjectives; i++) {
            for (int j = 0; j < numberObjectives - (i + 1); j++) {
                f[i] *= java.lang.Math.cos(java.lang.Math.pow(x[j], alpha) * (java.lang.Math.PI / 2.0));
            }
            if (i != 0) {
                int aux = numberObjectives - (i + 1);
                f[i] *= java.lang.Math.sin(java.lang.Math.pow(x[aux], alpha) * (java.lang.Math.PI / 2.0));
            }
        }

        return f;

    }
}
