package ec.app.TracableVectorProblems.MultiObjective.DTLZ;

import ec.EvolutionState;
import ec.app.TracableVectorProblems.MultiObjective.MultiObjectiveTraceableDoubleVectorProblem;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;

public class DTLZ7 extends MultiObjectiveTraceableDoubleVectorProblem {

    public void initialSetup(final EvolutionState state, final Parameter base) {
        this.minNumberObjectives = 2; //the number of objectives for this problem
        this.maxNumberObjectives = 3; //the number of objectives for this problem
        this.problemName = "DTLZ7"; //the name of the problem
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




        double g = 0.0;
        for (int i = genomeLength - k; i < genomeLength; i++) {
            g += x[i];
        }

        g = 1 + (9.0 * g) / k;

        System.arraycopy(x, 0, f, 0, numberObjectives - 1);

        double h = 0.0;
        for (int i = 0; i < numberObjectives - 1; i++) {
            h += (f[i] / (1.0 + g)) * (1 + Math.sin(3.0 * Math.PI * f[i]));
        }

        h = numberObjectives - h;

        f[numberObjectives - 1] = (1 + g) * h;

        return f;

    }
}
