package ec.vector;

import ec.util.MersenneTwisterFast;
import ec.EvolutionState;
import ec.Individual;
import ec.util.Code;
import ec.util.DecodeReturn;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableFloat;
import ec.vector.TracableDataTypes.TraceableInteger;
import ec.vector.TracableDataTypes.TraceTuple;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class TraceableIntegerVectorIndividual extends VectorIndividual {
    public static final String P_TRACABLEINTEGERINDIVIDUAL = "trace-bit-vect-ind";
    public TraceableInteger[] genome;

    private int randomID = -1;

    public Parameter defaultBase()
    {
        return VectorDefaults.base().push(P_TRACABLEINTEGERINDIVIDUAL);
    }

    public Object clone()
    {
        TraceableIntegerVectorIndividual myobj = (TraceableIntegerVectorIndividual) (super.clone());

        // must clone the genome
        //NOTE: needs to be this fancy to avoid reference errors
        myobj.genome = new TraceableInteger[genome.length];
        for(int i = 0; i < genome.length; i++){
            myobj.genome[i] = new TraceableInteger(genome[i].getValue(), genome[i].getTraceVector());
        }

        return myobj;
    }

    public void setup(final EvolutionState state, final Parameter base) //TODO: rework initialization!
    {
        super.setup(state,base);  // actually unnecessary (Individual.setup() is empty)

        IntegerVectorSpecies s = (IntegerVectorSpecies)species;  // where my default info is stored
        genome = new TraceableInteger[s.genomeSize];
        for (int i = 0; i < s.genomeSize; i++) {
            genome[i] = new TraceableInteger();
        }
    }

    /**
     * combines two TraceableFloat genes. Returns the new gene value.
     * @param a reference to the first gene
     * @param b reference to the second gene
     * @param new_a_value the new gene value of a
     * @return the new, combined gene
     */
    protected static TraceableInteger recombineGenes(TraceableInteger a, TraceableInteger b, int new_a_value){

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

        return new TraceableInteger(new_a_value, new_traceVector);
    }

    public void defaultCrossover(EvolutionState state, int thread, VectorIndividual ind)
    {
        IntegerVectorSpecies s = (IntegerVectorSpecies)species;  // where my default info is stored
        TraceableIntegerVectorIndividual i = (TraceableIntegerVectorIndividual) ind;
        TraceableInteger tmp;
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
            case VectorSpecies.C_TWO_POINT:
            {
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
                if (point0 > point) { int p = point0; point0 = point; point = p; }
                for(int x=point0*s.chunksize;x<point*s.chunksize;x++)
                {
                    tmp = i.genome[x];
                    i.genome[x] = genome[x];
                    genome[x] = tmp;
                }
            }
            break;
            case VectorSpecies.C_TWO_POINT_NO_NOP:
            {
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
            }
            break;
            case IntegerVectorSpecies.C_ANY_POINT:
                for(int x=0;x<len/s.chunksize;x++)
                    if (state.random[thread].nextBoolean(s.crossoverProbability))
                        for(int y=x*s.chunksize;y<(x+1)*s.chunksize;y++)
                        {
                            tmp = i.genome[y];
                            i.genome[y] = genome[y];
                            genome[y] = tmp;
                        }
                break;
            default:
                state.output.fatal("In valid crossover type in TracableIntegerVectorIndividual.");
                break;
            case VectorSpecies.C_LINE_RECOMB:
            {
                double alpha = state.random[thread].nextDouble() * (1 + 2*s.lineDistance) - s.lineDistance;
                double beta = state.random[thread].nextDouble() * (1 + 2*s.lineDistance) - s.lineDistance;
                long t,u;
                long min, max;
                for (int x = 0; x < len; x++)
                {
                    min = s.minGene(x);
                    max = s.maxGene(x);
                    t = (long) Math.floor(alpha * genome[x].getValue() + (1 - alpha) * i.genome[x].getValue() + 0.5);
                    u = (long) Math.floor(beta * i.genome[x].getValue() + (1 - beta) * genome[x].getValue() + 0.5);
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
                long t,u;
                long min, max;
                for (int x = 0; x < len; x++)
                {
                    do
                    {
                        double alpha = state.random[thread].nextDouble() * (1 + 2*s.lineDistance) - s.lineDistance;
                        double beta = state.random[thread].nextDouble() * (1 + 2*s.lineDistance) - s.lineDistance;
                        min = s.minGene(x);
                        max = s.maxGene(x);
                        t = (long) Math.floor(alpha * genome[x].getValue() + (1 - alpha) * i.genome[x].getValue() + 0.5);
                        u = (long) Math.floor(beta * i.genome[x].getValue() + (1 - beta) * genome[x].getValue() + 0.5);
                    } while (t < min || t > max || u < min || u > max);
                    genome[x] = this.recombineGenes(genome[x], i.genome[x], (int) t);
                    i.genome[x] = this.recombineGenes(i.genome[x], genome[x], (int) u);
                }
            }
            break;
        }
    }

    /** Splits the genome into n pieces, according to points, which *must* be sorted.
     pieces.length must be 1 + points.length */
    public void split(int[] points, Object[] pieces)
    {
        System.out.println("split called"); //DEBUGG: check if called, delete if i know when
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
        System.out.println("join called"); //DEBUGG: check if called, delete if i know when
        int sum=0;
        for(int x=0;x<pieces.length;x++)
            sum += ((boolean[])(pieces[x])).length;

        int runningsum = 0;
        TraceableInteger[] newgenome = new TraceableInteger[sum];
        for(int x=0;x<pieces.length;x++)
        {
            System.arraycopy(pieces[x], 0, newgenome, runningsum, ((boolean[])(pieces[x])).length);
            runningsum += ((boolean[])(pieces[x])).length;
        }
        // set genome
        genome = newgenome;
    }

    /**
     * sets the gene with 100% influence of the traceID
     * @param gene the gene to be set
     * @param value the value for the gene
     * @param traceID the traceID of the gene
     */
    private void setGene(TraceableInteger gene, int value, int traceID){
        List<TraceTuple> traceVector = new ArrayList<TraceTuple>();
        traceVector.add(new TraceTuple(traceID, 1.0));
        gene.setValue(value, traceVector);
    }

    /**
     * handles the mutation of one gene.
     * @param a the gene to be mutated
     * @param new_a_value the new value of gene a
     * @param mutationCounter the current mutation counter used for the traceID
     * @return the new, mutated gene
     */
    private TraceableInteger mutateGene(TraceableInteger a, int new_a_value, int mutationCounter){

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

        return new TraceableInteger(new_a_value, a_traceVector);
    }

    /** Destructively mutates the individual in some default manner.  The default form
     simply randomizes genes to a uniform distribution from the min and max of the gene values. */
    public void defaultMutate(EvolutionState state, int thread)
    {
        IntegerVectorSpecies s = (IntegerVectorSpecies) species;
        for(int x = 0; x < genome.length; x++)
            if (state.random[thread].nextBoolean(s.mutationProbability(x)))
            {
                int old = genome[x].getValue();
                for(int retries = 0; retries < s.duplicateRetries(x) + 1; retries++)
                {
                    switch(s.mutationType(x))
                    {
                        case IntegerVectorSpecies.C_RESET_MUTATION:
                        {
                            state.mutationCounter--;
                            //for random mutation the mutation value is completely random
                            int newValue = randomValueFromClosedInterval((int)s.minGene(x), (int)s.maxGene(x), state.random[thread]);
                            this.setGene(genome[x], newValue, state.mutationCounter);
                        }
                            break;
                        case IntegerVectorSpecies.C_RANDOM_WALK_MUTATION:
                            int min = (int)s.minGene(x);
                            int max = (int)s.maxGene(x);
                            int newValue = genome[x].getValue();
                            if (!s.mutationIsBounded(x))
                            {
                                // okay, technically these are still bounds, but we can't go beyond this without weird things happening
                                max = Integer.MAX_VALUE;
                                min = Integer.MIN_VALUE;
                            }
                            do
                            {
                                int n = (int)(state.random[thread].nextBoolean() ? 1 : -1);
                                if ((n == 1 && newValue < max) || (n == -1 && newValue > min)) {
                                    newValue = newValue + n;
                                }
                                else if ((n == -1 && newValue < max) || (n == 1 && newValue > min)) {
                                    newValue = newValue - n;
                                }
                            }
                            while (state.random[thread].nextBoolean(s.randomWalkProbability(x)));
                            state.mutationCounter--;
                            genome[x] = this.mutateGene(genome[x], newValue, state.mutationCounter);
                            break;
                        default:
                            state.output.fatal("In IntegerVectorIndividual.defaultMutate, default case occurred when it shouldn't have");
                            break;
                    }
                    if (genome[x].getValue() != old) break;
                    // else genome[x] = old;  // try again
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

    /** Initializes the individual by randomly assigning int values */
    public void reset(EvolutionState state, int thread)
    {
        IntegerVectorSpecies s = (IntegerVectorSpecies) species;
        for(int x=0;x<genome.length;x++) {
            List<TraceTuple> traceVector = new ArrayList<TraceTuple>();
            traceVector.add(new TraceTuple(this.traceID, 1.0));
            this.setGene(genome[x], randomValueFromClosedInterval((int) s.minGene(x), (int) s.maxGene(x), state.random[thread]), this.traceID);
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
        {
            if (i > 0) s.append("; ");
            s.append(genome[i].toString());
        }
        return s.toString();
    }

    /**
     * returns a string: [value, traceVector] => [value, [[traceID, impact,]...]
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

        genome = new TraceableInteger[ lll ];

        // read in the genes
        for( int i = 0 ; i < genome.length ; i++ )
        {
            Code.decode( d );
            genome[i] = new TraceableInteger();
            genome[i].fromString(d.s);
        }
    }

    public boolean equals(Object ind)
    {
        System.out.println("equals called"); //DEBUGG: check if called, delete if i know when
        if (ind==null) return false;
        if (!(this.getClass().equals(ind.getClass()))) return false; // SimpleRuleIndividuals are special.
        TraceableIntegerVectorIndividual i = (TraceableIntegerVectorIndividual)ind;
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
    { genome = (TraceableInteger[]) gen; }
    public int genomeLength()
    { return genome.length; }

    /** Clips each gene value to be within its specified [min,max] range. */
    public void clamp()
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
        System.out.println("setGenotypeLength called"); //DEBUGG: check if called, delete if i know when
        TraceableInteger[] newGenome = new TraceableInteger[len];
        System.arraycopy(genome, 0, newGenome, 0,
                genome.length < newGenome.length ? genome.length : newGenome.length);
        genome = newGenome;
    }

    public void writeGenotype(final EvolutionState state,
                              final DataOutput dataOutput) throws IOException
    {
        dataOutput.writeInt(genome.length);
        for(int x=0;x<genome.length;x++)
            dataOutput.writeBytes(genome[x].toString());  //TODO: fix!
            //dataOutput.writeBytes(genome[x].getValue() + "from individual "+ genome[x].getTraceVector().toString());  // is .writeBytes right? How is this method used?
    }

    public void readGenotype(final EvolutionState state,
                             final DataInput dataInput) throws IOException
    {
        System.out.println("read genotype called"); //DEBUGG: check if called, delete if i know when
            /*
            int len = dataInput.readInt();
            if (genome==null || genome.length != len)
                genome = new TracableInteger[len];
            for(int x=0;x<genome.length;x++)
                genome[x] = dataInput.readBoolean();
             */
        //TODO: implement
    }


    public double distanceTo(Individual otherInd)
    {
        TraceableIntegerVectorIndividual other = (TraceableIntegerVectorIndividual) otherInd;

        TraceableInteger[] otherGenome = other.genome;
        double sumSquaredDistance =0.0;
        for(int i=0; i < other.genomeLength(); i++)
        {
            long dist = this.genome[i].getValue() - (long)otherGenome[i].getValue();
            sumSquaredDistance += dist*dist;
        }
        return StrictMath.sqrt(sumSquaredDistance);
    }

}
