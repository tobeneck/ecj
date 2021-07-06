package ec.app.TracableVectorProblems.MultiObjective.DTLZ;

import ec.EvolutionState;
import ec.app.TracableVectorProblems.MultiObjective.MultiObjectiveTraceableDoubleVectorProblem;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;

public class DTLZ6 extends MultiObjectiveTraceableDoubleVectorProblem {

    public void initialSetup(final EvolutionState state, final Parameter base) {
        this.minNumberObjectives = 2; //the number of objectives for this problem
        this.maxNumberObjectives = 3; //the number of objectives for this problem
        this.problemName = "DTLZ6"; //the name of the problem
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



        double[] theta = new double[numberObjectives - 1];
        double g = 0.0;
        for (int i = genomeLength - k; i < genomeLength; i++) {
            g += java.lang.Math.pow(x[i], 0.1);
        }

        double t = java.lang.Math.PI / (4.0 * (1.0 + g));
        theta[0] = x[0] * java.lang.Math.PI / 2;
        for (int i = 1; i < (numberObjectives - 1); i++) {
            theta[i] = t * (1.0 + 2.0 * g * x[i]);
        }

        for (int i = 0; i < numberObjectives; i++) {
            f[i] = 1.0 + g;
        }

        for (int i = 0; i < numberObjectives; i++) {
            for (int j = 0; j < numberObjectives - (i + 1); j++) {
                f[i] *= java.lang.Math.cos(theta[j]);
            }
            if (i != 0) {
                int aux = numberObjectives - (i + 1);
                f[i] *= java.lang.Math.sin(theta[aux]);
            }
        }

        return f;

    }
}
