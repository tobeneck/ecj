package ec.vector;

import ec.util.MersenneTwisterFast;
import ec.EvolutionState;
import ec.Individual;
import ec.util.Code;
import ec.util.DecodeReturn;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceTuple;
import ec.vector.TracableDataTypes.TraceableDouble;

import java.util.ArrayList;
import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Arrays;
import java.util.List;

public class TraceableDoubleVectorIndividual extends VectorIndividual {

    private static final long serialVersionUID = 1;

    public static final String P_DOUBLEVECTORINDIVIDUAL = "double-vect-ind";

    public static final double MAXIMUM_INTEGER_IN_DOUBLE = 9.007199254740992E15;

    public static final String P_TRACEABLEDOUBLEVECTORINDIVIDUAL = "trace-double-vect-ind";
    public TraceableDouble[] genome;

    public Parameter defaultBase()
    {
        return VectorDefaults.base().push(P_TRACEABLEDOUBLEVECTORINDIVIDUAL);
    }

    public Object clone()
    {
        TraceableDoubleVectorIndividual myobj = (TraceableDoubleVectorIndividual) (super.clone());

        // must clone the genome
        //NOTE: needs to be this fancy to avoid reference errors
        myobj.genome = new TraceableDouble[genome.length];
        for(int i = 0; i < genome.length; i++){
            myobj.genome[i] = new TraceableDouble(genome[i].getValue(), genome[i].getTraceVector());
        }

        return myobj;
    }

    public void setup(final EvolutionState state, final Parameter base)
    {
        super.setup(state,base);  // actually unnecessary (Individual.setup() is empty)

        FloatVectorSpecies s = (FloatVectorSpecies)species;  // where my default info is stored
        genome = new TraceableDouble[s.genomeSize];
        for (int i = 0; i < s.genomeSize; i++) {
            genome[i] = new TraceableDouble();
        }
    }

    //TODO: does this inline stuff work?
    public void simulatedBinaryCrossover(MersenneTwisterFast random, TraceableDoubleVectorIndividual other, double eta_c)
    {
        final double EPS = FloatVectorSpecies.SIMULATED_BINARY_CROSSOVER_EPS;
        FloatVectorSpecies s = (FloatVectorSpecies) species;
        TraceableDouble[] parent1 = genome;
        TraceableDouble[] parent2 = other.genome;
        //        double[] min_realvar = s.minGenes;
        //        double[] max_realvar = s.maxGenes;


        double y1, y2, yl, yu;
        double c1, c2;
        double alpha, beta, betaq;
        double rand;

        for(int i = 0; i < parent1.length; i++)
        {
            if (random.nextBoolean())  // 0.5f
            {
                if (Math.abs(parent1[i].getValue() - parent2[i].getValue()) > EPS)
                {
                    if (parent1[i].getValue() < parent2[i].getValue())
                    {
                        y1 = parent1[i].getValue();
                        y2 = parent2[i].getValue();
                    }
                    else
                    {
                        y1 = parent2[i].getValue();
                        y2 = parent1[i].getValue();
                    }
                    yl = s.minGene(i); //min_realvar[i];
                    yu = s.maxGene(i); //max_realvar[i];
                    rand = random.nextDouble();
                    beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                    alpha = 2.0 - Math.pow(beta,-(eta_c+1.0));
                    if (rand <= (1.0/alpha))
                    {
                        betaq = Math.pow((rand*alpha),(1.0/(eta_c+1.0)));
                    }
                    else
                    {
                        betaq = Math.pow((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                    }
                    c1 = 0.5*((y1+y2)-betaq*(y2-y1));
                    beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                    alpha = 2.0 - Math.pow(beta,-(eta_c+1.0));
                    if (rand <= (1.0/alpha))
                    {
                        betaq = Math.pow((rand*alpha),(1.0/(eta_c+1.0)));
                    }
                    else
                    {
                        betaq = Math.pow((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                    }
                    c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                    if (c1<yl)
                        c1=yl;
                    if (c2<yl)
                        c2=yl;
                    if (c1>yu)
                        c1=yu;
                    if (c2>yu)
                        c2=yu;
                    if (random.nextBoolean())
                    {
                        parent1[i] = this.recombineGenes(parent1[i], parent2[i], (double) c2);
                        parent2[i] = this.recombineGenes(parent2[i], parent1[i], (double) c1);
                    }
                    else
                    {
                        parent1[i] = this.recombineGenes(parent1[i], parent2[i], (double) c2);
                        parent2[i] = this.recombineGenes(parent2[i], parent1[i], (double) c1);
                    }
                }
                else
                {
                    // do nothing
                }
            }
            else
            {
                // do nothing
            }
        }
    }


    /**
     * combines two TraceableDouble genes. Returns the new gene value.
     * @param a reference to the first gene
     * @param b reference to the second gene
     * @param new_a_value the new gene value of a
     * @return the new, combined gene
     */
    protected static TraceableDouble recombineGenes(TraceableDouble a, TraceableDouble b, double new_a_value){

        // of aVal = bVal = newVal, both genes should have 50% impact
        double influence_factor_a = 0.5;
        double influence_factor_b = 0.5;

        if(a.getValue() == new_a_value && b.getValue() != new_a_value) //return a if it has 100% influence. Faster and avoids 0.0 impact traceTuples
            return a;
        if(b.getValue() == new_a_value && a.getValue() != new_a_value) //return b if it has 100% influence. Faster and avoids 0.0 impact traceTuples
            return b;
        if(a.getValue() != new_a_value && b.getValue() != new_a_value) { //compute the new values if non of them are equal
            //if you don't cast here, the resulting values are integers and will be roundet!
            influence_factor_a = 1.0 - (double)Math.abs(a.getValue() - new_a_value) / (double)(Math.abs(a.getValue() - new_a_value) + Math.abs(b.getValue() - new_a_value));
            influence_factor_b = 1.0 - (double)Math.abs(b.getValue() - new_a_value) / (double)(Math.abs(a.getValue() - new_a_value) + Math.abs(b.getValue() - new_a_value));
        }

        int i = 0; //index for trace vector a
        int j = 0; //index for trace vector b

        List<TraceTuple> new_traceVector = new ArrayList<TraceTuple>();

        while(true){ //this iterates over the traceVector of this individual
            if(i >= a.getTraceVector().size() && j >= b.getTraceVector().size()) //stop if both vectors are empty
                break;
            else if(i >= a.getTraceVector().size() && !(j >= b.getTraceVector().size())){//append if the a vector is empty and b vector is not.
                int currentBID = b.getTraceVector().get(j).getTraceID();
                double currentBImpact = b.getTraceVector().get(j).getImpact();
                new_traceVector.add(new TraceTuple(currentBID, influence_factor_b * currentBImpact));
                j++;
            }
            else if(!(i >= a.getTraceVector().size()) && j >= b.getTraceVector().size()){//append if the b vector is empty and a vector is not.
                int currentAID = a.getTraceVector().get(i).getTraceID();
                double currentAImpact = a.getTraceVector().get(i).getImpact();
                new_traceVector.add(new TraceTuple(currentAID, influence_factor_a * currentAImpact));
                i++;
            }
            else {//if both arrays are not empty, append the next traceID:
                int currentAID = a.getTraceVector().get(i).getTraceID();
                int currentBID = b.getTraceVector().get(j).getTraceID();

                double currentAImpact = a.getTraceVector().get(i).getImpact();
                double currentBImpact = b.getTraceVector().get(j).getImpact();

                if (currentAID == currentBID) {//combine the two if equal
                    new_traceVector.add(new TraceTuple(currentAID, influence_factor_a * currentAImpact + influence_factor_b * currentBImpact));
                    i++;
                    j++;
                }

                if (currentAID < currentBID) {//add the traceID of a if its smaller than the traceID of b
                    new_traceVector.add(new TraceTuple(currentAID, influence_factor_a * currentAImpact));
                    i++;
                }

                if (currentBID < currentAID) {//add the traceID of b if its smaller than the traceID of a
                    new_traceVector.add(new TraceTuple(currentBID, influence_factor_b * currentBImpact));
                    j++;
                }
            }
        }

        return new TraceableDouble(new_a_value, new_traceVector);
    }

    public void defaultCrossover(EvolutionState state, int thread, VectorIndividual ind)
    {
        FloatVectorSpecies s = (FloatVectorSpecies)species;  // where my default info is stored
        TraceableDoubleVectorIndividual i = (TraceableDoubleVectorIndividual) ind;
        TraceableDouble tmp;
        int point;

        int len = Math.min(genome.length, i.genome.length);
        if (len != genome.length || len != i.genome.length)
            state.output.warnOnce("Genome lengths are not the same.  Vector crossover will only be done in overlapping region.");

        switch(s.crossoverType)
        {
            case VectorSpecies.C_ONE_POINT:
                //                point = state.random[thread].nextInt((len / s.chunksize)+1);
                // we want to go from 0 ... len-1
                // so that there is only ONE case of NO-OP crossover, not TWO
                point = state.random[thread].nextInt((len / s.chunksize));
                for(int x=0;x<point*s.chunksize;x++)
                {
                    tmp = i.genome[x];
                    i.genome[x] = genome[x];
                    genome[x] = tmp;
                }
                break;
            case VectorSpecies.C_ONE_POINT_NO_NOP:
                point = state.random[thread].nextInt((len / s.chunksize) - 1) + 1;  // so it goes from 1 .. len-1
                for(int x=0;x<point*s.chunksize;x++)
                {
                    tmp = i.genome[x];
                    i.genome[x] = genome[x];
                    genome[x] = tmp;
                }
                break;
            case VectorSpecies.C_TWO_POINT: {
                //                int point0 = state.random[thread].nextInt((len / s.chunksize)+1);
                //                point = state.random[thread].nextInt((len / s.chunksize)+1);
                // we want to go from 0 to len-1
                // so that the only NO-OP crossover possible is point == point0
                // example; len = 4
                // possibilities: a=0 b=0       NOP                             [0123]
                //                                a=0 b=1       swap 0                  [for 1, 2, 3]
                //                                a=0 b=2       swap 0, 1               [for 2, 3]
                //                                a=0 b=3       swap 0, 1, 2    [for 3]
                //                                a=1 b=1       NOP                             [1230]
                //                                a=1 b=2       swap 1                  [for 2, 3, 0]
                //                                a=1 b=3       swap 1, 2               [for 3, 0]
                //                                a=2 b=2       NOP                             [2301]
                //                                a=2 b=3       swap 2                  [for 3, 0, 1]
                //                                a=3 b=3   NOP                         [3012]
                // All intervals: 0, 01, 012, 0123, 1, 12, 123, 1230, 2, 23, 230, 2301, 3, 30, 301, 3012
                point = state.random[thread].nextInt((len / s.chunksize));
                int point0 = state.random[thread].nextInt((len / s.chunksize));
                if (point0 > point) {
                    int p = point0;
                    point0 = point;
                    point = p;
                }
                for (int x = point0 * s.chunksize; x < point * s.chunksize; x++) {
                    tmp = i.genome[x];
                    i.genome[x] = genome[x];
                    genome[x] = tmp;
                }
            }
            break;
            case VectorSpecies.C_TWO_POINT_NO_NOP:
                point = state.random[thread].nextInt((len / s.chunksize));
                int point0 = 0;
                do { point0 = state.random[thread].nextInt((len / s.chunksize)); }
                while (point0 == point);  // NOP
                if (point0 > point) { int p = point0; point0 = point; point = p; }
                for(int x=point0*s.chunksize;x<point*s.chunksize;x++)
                {
                    tmp = i.genome[x];
                    i.genome[x] = genome[x];
                    genome[x] = tmp;
                }
                break;
            case VectorSpecies.C_ANY_POINT:
                for(int x=0;x<len/s.chunksize;x++)
                    if (state.random[thread].nextBoolean(s.crossoverProbability))
                        for(int y=x*s.chunksize;y<(x+1)*s.chunksize;y++)
                        {
                            tmp = i.genome[y];
                            i.genome[y] = genome[y];
                            genome[y] = tmp;
                        }
                break;
            case VectorSpecies.C_LINE_RECOMB:
            {
                double alpha = state.random[thread].nextDouble(true, true) * (1 + 2*s.lineDistance) - s.lineDistance;
                double beta = state.random[thread].nextDouble(true, true) * (1 + 2*s.lineDistance) - s.lineDistance;
                double t,u,min,max;
                for (int x = 0; x < len; x++)
                {
                    min = s.minGene(x);
                    max = s.maxGene(x);
                    t = alpha * genome[x].getValue() + (1 - alpha) * i.genome[x].getValue();
                    u = beta * i.genome[x].getValue() + (1 - beta) * genome[x].getValue();
                    if (!(t < min || t > max || u < min || u > max))
                    {
                        genome[x] = this.recombineGenes(genome[x], i.genome[x], (int) t);
                        i.genome[x] = this.recombineGenes(i.genome[x], genome[x], (int) u);
                    }
                }
            }
            break;
            case VectorSpecies.C_INTERMED_RECOMB:
            {
                double t,u,min,max;
                for (int x = 0; x < len; x++)
                {
                    do
                    {
                        double alpha = state.random[thread].nextDouble(true, true) * (1 + 2*s.lineDistance) - s.lineDistance;
                        double beta = state.random[thread].nextDouble(true, true) * (1 + 2*s.lineDistance) - s.lineDistance;
                        min = s.minGene(x);
                        max = s.maxGene(x);
                        t = alpha * genome[x].getValue() + (1 - alpha) * i.genome[x].getValue();
                        u = beta * i.genome[x].getValue() + (1 - beta) * genome[x].getValue();
                    } while (t < min || t > max || u < min || u > max);
                    genome[x] = this.recombineGenes(genome[x], i.genome[x], (int) t);
                    i.genome[x] = this.recombineGenes(i.genome[x], genome[x], (int) u);
                }
            }
            break;
            case VectorSpecies.C_SIMULATED_BINARY:
                simulatedBinaryCrossover(state.random[thread], i, s.crossoverDistributionIndex);
                break;
            default:
                state.output.fatal("In valid crossover type in TraceableDoubleVectorIndividual.");
                break;
        }
    }

    /** Splits the genome into n pieces, according to points, which *must* be sorted.
     pieces.length must be 1 + points.length */
    public void split(int[] points, Object[] pieces)
    {
        int point0, point1;
        point0 = 0; point1 = points[0];
        for(int x=0;x<pieces.length;x++)
        {
            pieces[x] = new boolean[point1-point0];
            System.arraycopy(genome,point0,pieces[x],0,point1-point0);
            point0 = point1;
            if (x >=pieces.length-2)
                point1 = genome.length;
            else point1 = points[x+1];
        }
    }

    /** Joins the n pieces and sets the genome to their concatenation.*/
    public void join(Object[] pieces)
    {
        int sum=0;
        for(int x=0;x<pieces.length;x++)
            sum += ((boolean[])(pieces[x])).length;

        int runningsum = 0;
        TraceableDouble[] newgenome = new TraceableDouble[sum];
        for(int x=0;x<pieces.length;x++)
        {
            System.arraycopy(pieces[x], 0, newgenome, runningsum, ((boolean[])(pieces[x])).length);
            runningsum += ((boolean[])(pieces[x])).length;
        }
        // set genome
        genome = newgenome;
    }

    /**
     * handles the mutation of one gene.
     * @param a the gene to be mutated
     * @param new_a_value the new value of gene a
     * @param mutationCounter the current mutation counter used for the traceID
     * @return the new, mutated gene
     */
    private TraceableDouble mutateGene(TraceableDouble a, double new_a_value, int mutationCounter){

        if(a.getValue() == new_a_value)//return the original gene if the value was not altered
            return a;

        double influence_factor_mut = Math.abs(a.getValue() - new_a_value) / (Math.abs(a.getValue()) + Math.abs(a.getValue() - new_a_value));
        double influence_factor_old = 1 - influence_factor_mut;

        List<TraceTuple> a_traceVector = new ArrayList<TraceTuple>();

        a_traceVector.add(new TraceTuple(mutationCounter, influence_factor_mut));

        int i = 0; //index for this traceVector
        while(i < a.getTraceVector().size()){ //this iterates over the traceVector of this individual
            int currentAID = a.getTraceVector().get(i).getTraceID();
            int currentAImpact = a.getTraceVector().get(i).getTraceID();

            a_traceVector.add(new TraceTuple(currentAID, influence_factor_old * currentAImpact));
            i++;
        }

        return new TraceableDouble(new_a_value, a_traceVector);
    }

    /**
     * Destructively mutates the individual in some default manner. The default
     * form simply randomizes genes to a uniform distribution from the min and
     * max of the gene values. It can also add gaussian noise to the genes, if
     * so directed in the FloatVectorSpecies. If the gaussian noise pushes the
     * gene out of range, a new noise value is generated.
     *
     * @author Sean Luke, Liviu Panait and Gabriel Balan
     */
    public void defaultMutate(EvolutionState state, int thread)
    {
        FloatVectorSpecies s = (FloatVectorSpecies) species;

        MersenneTwisterFast rng = state.random[thread];
        for(int x = 0; x < genome.length; x++)
            if (rng.nextBoolean(s.mutationProbability(x)))
            {
                double old = genome[x].getValue();
                for(int retries = 0; retries < s.duplicateRetries(x) + 1; retries++)
                {
                    switch(s.mutationType(x))
                    {
                        case FloatVectorSpecies.C_GAUSS_MUTATION:
                            gaussianMutation(state, rng, s, x);
                            break;
                        case FloatVectorSpecies.C_POLYNOMIAL_MUTATION:
                            polynomialMutation(state, rng, s, x);
                            break;
                        case FloatVectorSpecies.C_RESET_MUTATION:
                            floatResetMutation(state, rng, s, x);
                            break;
                        case FloatVectorSpecies.C_INTEGER_RESET_MUTATION:
                            integerResetMutation(state, rng, s, x);
                            break;
                        case FloatVectorSpecies.C_INTEGER_RANDOM_WALK_MUTATION:
                            integerRandomWalkMutation(state, rng, s, x);
                            break;
                        default:
                            state.output.fatal("In TraceableDoubleVectorIndividual.defaultMutate, default case occurred when it shouldn't have");
                            break;
                    }
                    if (genome[x].getValue() != old) break;
                    // else genome[x] = old;  // try again
                }
            }
    }

    void integerRandomWalkMutation(EvolutionState state, MersenneTwisterFast random, FloatVectorSpecies species, int index)
    {
        double min = species.minGene(index);
        double max = species.maxGene(index);
        if (!species.mutationIsBounded(index))
        {
            // okay, technically these are still bounds, but we can't go beyond this without weird things happening
            max = MAXIMUM_INTEGER_IN_DOUBLE;
            min = -(max);
        }
        double newValue = (double)Math.floor(genome[index].getValue());
        do
        {
            int n = (int)(random.nextBoolean() ? 1 : -1);
            if ((n == 1 && newValue < max) ||
                    (n == -1 && newValue > min)) {
                newValue = n + newValue;
            }
            else if ((n == -1 && newValue < max) ||
                    (n == 1 && newValue > min)) {
                newValue = n - newValue;
            }
        }
        while (random.nextBoolean(species.randomWalkProbability(index)));

        state.mutationCounter--;
        genome[index] = this.mutateGene(genome[index], newValue, state.mutationCounter);
    }

    void integerResetMutation(EvolutionState state, MersenneTwisterFast random, FloatVectorSpecies species, int index)
    {
        int minGene = (int)Math.floor(species.minGene(index));
        int maxGene = (int)Math.floor(species.maxGene(index));
        state.mutationCounter--;
        genome[index] = new TraceableDouble(randomValueFromClosedInterval(minGene, maxGene, random), state.mutationCounter);// minGene + random.nextLong(maxGene - minGene + 1);
    }

    void floatResetMutation(EvolutionState state, MersenneTwisterFast random, FloatVectorSpecies species, int index)
    {
        double minGene = species.minGene(index);
        double maxGene = species.maxGene(index);
        state.mutationCounter--;
        genome[index] = new TraceableDouble((double)(minGene + random.nextDouble(true, true) * (maxGene - minGene)), state.mutationCounter);
    }

    void gaussianMutation(EvolutionState state, MersenneTwisterFast random, FloatVectorSpecies species, int index)
    {
        double val;
        double min = species.minGene(index);
        double max = species.maxGene(index);
        double stdev = species.gaussMutationStdev(index);
        int outOfBoundsLeftOverTries = species.outOfBoundsRetries;
        boolean givingUpAllowed = species.outOfBoundsRetries != 0;
        do
        {
            val = random.nextGaussian() * stdev + genome[index].getValue();
            outOfBoundsLeftOverTries--;
            if (species.mutationIsBounded(index) && (val > max || val < min))
            {
                if (givingUpAllowed && (outOfBoundsLeftOverTries == 0))
                {
                    val = min + random.nextDouble() * (max - min);
                    species.outOfRangeRetryLimitReached(state);// it better get inlined
                    break;
                }
            }
            else break;
        }
        while (true);

        state.mutationCounter--;
        genome[index] = this.mutateGene(genome[index], (double)val, state.mutationCounter);
    }

    void polynomialMutation(EvolutionState state, MersenneTwisterFast random, FloatVectorSpecies species, int index)
    {
        double eta_m = species.mutationDistributionIndex(index);
        boolean alternativePolynomialVersion = species.polynomialIsAlternative(index);

        double rnd, delta1, delta2, mut_pow, deltaq;
        double y, yl, yu, val, xy;
        double y1;

        y1 = y = genome[index].getValue();  // ind[index];
        yl = species.minGene(index); // min_realvar[index];
        yu = species.maxGene(index); // max_realvar[index];
        delta1 = (y-yl)/(yu-yl);
        delta2 = (yu-y)/(yu-yl);

        int totalTries = species.outOfBoundsRetries;
        int tries = 0;
        for(tries = 0; tries < totalTries || totalTries == 0; tries++)  // keep trying until totalTries is reached if it's not zero.  If it's zero, go on forever.
        {
            rnd = random.nextDouble();
            mut_pow = 1.0/(eta_m+1.0);
            if (rnd <= 0.5)
            {
                xy = 1.0-delta1;
                val = 2.0*rnd + (alternativePolynomialVersion ? (1.0-2.0*rnd)*(Math.pow(xy,(eta_m+1.0))) : 0.0);
                deltaq =  Math.pow(val,mut_pow) - 1.0;
            }
            else
            {
                xy = 1.0-delta2;
                val = 2.0*(1.0-rnd) + (alternativePolynomialVersion ? 2.0*(rnd-0.5)*(Math.pow(xy,(eta_m+1.0))) : 0.0);
                deltaq = 1.0 - (Math.pow(val,mut_pow));
            }
            y1 = y + deltaq*(yu-yl);
            if (!species.mutationIsBounded(index) || (y1 >= yl && y1 <= yu)) break;  // yay, found one
        }

        // at this point, if tries is totalTries, we failed
        if (totalTries != 0 && tries == totalTries)
        {
            // just randomize
            y1 = (double)(species.minGene(index) + random.nextDouble(true, true) * (species.maxGene(index) - species.minGene(index)));  //(double)(min_realvar[index] + random.nextDouble
            // () * (max_realvar[index] - min_realvar[index]));
            species.outOfRangeRetryLimitReached(state);// it better get inlined
        }

        state.mutationCounter--;
        genome[index] = this.mutateGene(genome[index], (double)y1, state.mutationCounter); // ind[index] = y1;
    }

    /** This function is broken out to keep it identical to NSGA-II's mutation.c code. eta_m is the distribution
     index.  */
    public void polynomialMutate(EvolutionState state, MersenneTwisterFast random, double eta_m, boolean alternativePolynomialVersion, boolean mutationIsBounded)
    {
        FloatVectorSpecies s = (FloatVectorSpecies) species;
        TraceableDouble[] ind = genome;
        //        double[] min_realvar = s.minGenes;
        //        double[] max_realvar = s.maxGenes;

        double rnd, delta1, delta2, mut_pow, deltaq;
        double y, yl, yu, val, xy;
        double y1;
        for (int j=0; j < ind.length; j++)
        {
            if (random.nextBoolean(s.mutationProbability[j]))
            {
                y1 = y = ind[j].getValue();
                yl = s.minGene(j); //min_realvar[j];
                yu = s.maxGene(j); //max_realvar[j];
                delta1 = (y-yl)/(yu-yl);
                delta2 = (yu-y)/(yu-yl);

                int totalTries = s.outOfBoundsRetries;
                int tries = 0;
                for(tries = 0; tries < totalTries || totalTries == 0; tries++)  // keep trying until totalTries is reached if it's not zero.  If it's zero, go on forever.
                {
                    rnd = (random.nextDouble());
                    mut_pow = 1.0/(eta_m+1.0);
                    if (rnd <= 0.5)
                    {
                        xy = 1.0-delta1;
                        val = 2.0*rnd + (alternativePolynomialVersion ? (1.0-2.0*rnd)*(Math.pow(xy,(eta_m+1.0))) : 0.0);
                        deltaq =  Math.pow(val,mut_pow) - 1.0;
                    }
                    else
                    {
                        xy = 1.0-delta2;
                        val = 2.0*(1.0-rnd) + (alternativePolynomialVersion ? 2.0*(rnd-0.5)*(Math.pow(xy,(eta_m+1.0))) : 0.0);
                        deltaq = 1.0 - (Math.pow(val,mut_pow));
                    }
                    y1 = y + deltaq*(yu-yl);
                    if (!mutationIsBounded || (y1 >= yl && y1 <= yu)) break;  // yay, found one
                }

                // at this point, if tries is totalTries, we failed
                if (totalTries != 0 && tries == totalTries)
                {
                    // just randomize
                    // y1 = (double)(min_realvar[j] + random.nextDouble(true, true) * (max_realvar[j] - min_realvar[j]));
                    y1 = (double)(s.minGene(j) + random.nextDouble(true, true) * (s.maxGene(j) - s.minGene(j)));
                    s.outOfRangeRetryLimitReached(state);// it better get inlined
                }

                state.mutationCounter--;
                ind[j] = this.mutateGene(ind[j], (double)y1, state.mutationCounter);
            }
        }
    }


    /** Returns a random value from between min and max inclusive.  This method handles
     overflows that complicate this computation.  Does NOT check that
     min is less than or equal to max.  You must check this yourself. */
    public int randomValueFromClosedInterval(int min, int max, MersenneTwisterFast random)
    {
        if (max - min < 0) // we had an overflow
        {
            int l = 0;
            do l = random.nextInt();
            while(l < min || l > max);
            return l;
        }
        else return min + random.nextInt(max - min + 1);
    }

    /**
     * Initializes the individual by randomly choosing doubles uniformly from
     * mingene to maxgene.
     */
    public void reset(EvolutionState state, int thread)
    {
        FloatVectorSpecies s = (FloatVectorSpecies) species;
        MersenneTwisterFast random = state.random[thread];
        for (int x = 0; x < genome.length; x++)
        {
            int type = s.mutationType(x);
            if (type == FloatVectorSpecies.C_INTEGER_RESET_MUTATION ||
                    type == FloatVectorSpecies.C_INTEGER_RANDOM_WALK_MUTATION)  // integer type
            {
                int minGene = (int)Math.floor(s.minGene(x));
                int maxGene = (int)Math.floor(s.maxGene(x));
                genome[x] = new TraceableDouble(randomValueFromClosedInterval(minGene, maxGene, random), this.traceID);//minGene + random.nextInt(maxGene - minGene + 1);
            }
            else
            {
                genome[x] = new TraceableDouble((double)(s.minGene(x) + random.nextDouble(true, true) * (s.maxGene(x) - s.minGene(x))), this.traceID);
            }
        }
    }

    public int hashCode()
    {
        // stolen from GPIndividual.  It's a decent algorithm.
        int hash = this.getClass().hashCode();

        hash = ( hash << 1 | hash >>> 31 ) ^ Arrays.hashCode(genome);

        return hash;
    }

    public String genotypeToStringForHumans()
    {
        StringBuilder s = new StringBuilder();
        for( int i = 0 ; i < genome.length ; i++ )
        { if (i > 0) s.append("; "); s.append(genome[i].toString()); }
        return s.toString();
    }

    /**
     * returns a string: value; traceID
     * @return
     */
    public String genotypeToString()
    {
        StringBuilder s = new StringBuilder();
        s.append( Code.encode( genome.length ) );
        for( int i = 0 ; i < genome.length ; i++ ) {
            s.append(Code.encode(genome[i].toString()));
        }
        return s.toString();
    }

    protected void parseGenotype(final EvolutionState state,
                                 final LineNumberReader reader) throws IOException
    {
        String s = reader.readLine();
        DecodeReturn d = new DecodeReturn(s);
        Code.decode( d );

        int lll = (int)(d.l);

        genome = new TraceableDouble[ lll ];

        // read in the genes
        for( int i = 0 ; i < genome.length ; i++ )
        {
            Code.decode( d );
            genome[i] = new TraceableDouble();
            genome[i].fromString(d.s);
        }
    }

    public boolean equals(Object ind)
    {
        System.out.println("equals called"); //DEBUGG: check if called, delete if i know when
        if (ind==null) return false;
        if (!(this.getClass().equals(ind.getClass()))) return false; // SimpleRuleIndividuals are special.
        TraceableDoubleVectorIndividual i = (TraceableDoubleVectorIndividual)ind;
        if( genome.length != i.genome.length )
            return false;
        for( int j = 0 ; j < genome.length ; j++ )
            if( genome[j] != i.genome[j] )
                return false;
        return true;
    }

    public Object getGenome()
    { return genome; }
    public void setGenome(Object gen)
    { genome = (TraceableDouble[]) gen; }
    public int genomeLength()
    { return genome.length; }

    /** Clips each gene value to be within its specified [min,max] range. */
    public void clamp() //TODO: implement
    {
        System.out.println("clamp was called!"); //DEBUGG: check if called, delete if i know when
        /*
        IntegerVectorSpecies _species = (IntegerVectorSpecies)species;
        for (int i = 0; i < genomeLength(); i++)
        {
            int minGene = (int)_species.minGene(i);
            if (genome[i].getValue() < minGene)
                genome[i].setValue(minGene); //what do do for the TraceID?
            else
            {
                int maxGene = (int)_species.maxGene(i);
                if (genome[i].getValue() > maxGene)
                    genome[i] = maxGene;
            }
        }

         */
    }

    public void setGenomeLength(int len)
    {
        TraceableDouble[] newGenome = new TraceableDouble[len];
        System.arraycopy(genome, 0, newGenome, 0,
                genome.length < newGenome.length ? genome.length : newGenome.length);
        genome = newGenome;
    }

    public void writeGenotype(final EvolutionState state,
                              final DataOutput dataOutput) throws IOException
    {
        dataOutput.writeInt(genome.length);
        for(int x=0;x<genome.length;x++)
            dataOutput.writeBytes(genome[x].getValue() + "from individual "+ genome[x].getTraceVector().toString());  // is .writeBytes right? How is this method used?
    }

    public void readGenotype(final EvolutionState state,
                             final DataInput dataInput) throws IOException
    {
        System.out.println("read genotype called"); //DEBUGG: check if called, delete if i know when
            /*
            int len = dataInput.readInt();
            if (genome==null || genome.length != len)
                genome = new TraceableDouble[len];
            for(int x=0;x<genome.length;x++)
                genome[x] = dataInput.readBoolean();
             */
    }


    public double distanceTo(Individual otherInd)
    {
        TraceableDoubleVectorIndividual other = (TraceableDoubleVectorIndividual) otherInd;

        TraceableDouble[] otherGenome = other.genome;
        double sumSquaredDistance =0.0;
        for(int i=0; i < other.genomeLength(); i++)
        {
            double dist = this.genome[i].getValue() - (long)otherGenome[i].getValue();
            sumSquaredDistance += dist*dist;
        }
        return StrictMath.sqrt(sumSquaredDistance);
    }

}
