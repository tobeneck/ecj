package ec.app.TracableVectorProblems.MultiObjective.ZDT;

import ec.EvolutionState;
import ec.app.TracableVectorProblems.MultiObjective.MultiObjectiveTraceableDoubleVectorProblem;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableDouble;

public class ZDT1 extends MultiObjectiveTraceableDoubleVectorProblem {

    public void initialSetup(final EvolutionState state, final Parameter base) {
        this.minNumberObjectives = 2; //the number of objectives for this problem
        this.maxNumberObjectives = 2; //the number of objectives for this problem
        this.problemName = "ZDT1"; //the name of the problem
    }

    /**
     * Returns the value of the ZDT1 function G.
     *
     * @param genome Genome
     */
    protected double evalG(TraceableDouble[] genome) {
        double g = 0.0;
        for (int i = 1; i < genome.length; i++) {
            g += genome[i].getValue();
        }
        double constant = 9.0 / (genome.length - 1); //TODO: shouldn't this be 29?

        return constant * g + 1.0;
    }

    /**
     * Returns the value of the ZDT1 function H.
     *
     * @param f1 First argument of the function H.
     * @param g Second argument of the function H.
     */
    protected double evalH(double f1, double g) {
        double h ;
        h = 1.0 - Math.sqrt(f1 / g);
        return h;
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

        f[0] = genome[0].getValue();
        double g = this.evalG(genome);
        double h = this.evalH(f[0], g);
        f[1] = h * g;

        return f;

    }

}
