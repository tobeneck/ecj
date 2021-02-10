package ec.vector;

import ec.vector.TracableDataTypes.TraceTuple;
import ec.vector.TracableDataTypes.TraceableBoolean;
import ec.EvolutionState;
import ec.Individual;
import ec.util.Code;
import ec.util.DecodeReturn;
import ec.util.Parameter;
import ec.vector.TracableDataTypes.TraceableFloat;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class TraceableBitVectorIndividual extends VectorIndividual {
    public static final String P_TRACABLEBITVECTORINDIVIDUAL = "trace-bit-vect-ind";
    public TraceableBoolean[] genome;

    private int randomID = -1;

    public Parameter defaultBase()
    {
        System.out.println("defaultBase called"); //DEBUGG: check if called, delete if i know when
        return VectorDefaults.base().push(P_TRACABLEBITVECTORINDIVIDUAL);
    }

    public Object clone()
    {
        TraceableBitVectorIndividual myobj = (TraceableBitVectorIndividual) (super.clone());

        // must clone the genome
        //NOTE: needs to be this fancy to avoid reference errors
        myobj.genome = new TraceableBoolean[genome.length];
        for(int i = 0; i < genome.length; i++){
            myobj.genome[i] = new TraceableBoolean(genome[i].getValue(), genome[i].getTraceVector());
        }

        return myobj;
    }

    public void setup(final EvolutionState state, final Parameter base) //TODO: rework initialization!
    {
        super.setup(state,base);  // actually unnecessary (Individual.setup() is empty)

        Parameter def = defaultBase();

        if (!(species instanceof BitVectorSpecies))
            state.output.fatal("TracableBitVectorIndividual requires an BitVectorSpecies", base, def);
        BitVectorSpecies s = (BitVectorSpecies) species;

        genome = new TraceableBoolean[s.genomeSize];
        for (int i = 0; i < s.genomeSize; i++) {
            genome[i] = new TraceableBoolean();
        }
    }

    public void defaultCrossover(EvolutionState state, int thread, VectorIndividual ind)
    {
        BitVectorSpecies s = (BitVectorSpecies)species;  // where my default info is stored
        TraceableBitVectorIndividual i = (TraceableBitVectorIndividual) ind;
        TraceableBoolean tmp;
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
                {//means, the gene

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
            case BitVectorSpecies.C_ANY_POINT:
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
                state.output.fatal("In valid crossover type in TracableBitVectorIndividual.");
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
        TraceableBoolean[] newgenome = new TraceableBoolean[sum];
        for(int x=0;x<pieces.length;x++)
        {
            System.arraycopy(pieces[x], 0, newgenome, runningsum, ((boolean[])(pieces[x])).length);
            runningsum += ((boolean[])(pieces[x])).length;
        }
        // set genome
        genome = newgenome;
    }

    /**
     * handles the mutation of one gene. Sets both the value and the traceVector.
     * Disclaimer: unlike the other traceable individuals, this method only works with bit flip mutation. The old and new gene values are always combined 50/50.
     * @param a the gene to be mutated
     * @param new_a_value the new value of gene a
     */
    private void mutateGene(TraceableBoolean a, boolean new_a_value, int mutationCounter){
        //as the value is 100% randomly generated, there is no need to compare the old gene values.

        double influence_factor_mut = 0.5; //as there are just two values, the influence of a mutation is always 50%
        double influence_factor_old = 1 - influence_factor_mut;

        List<TraceTuple> a_traceVector = new ArrayList<TraceTuple>();

        a_traceVector.add(new TraceTuple(mutationCounter, influence_factor_mut));

        int i = 0; //index for this traceVector
        while(i < a.getTraceVector().size()){ //this iterates over the traceVector of this individual
            int currentAID = a.getTraceVector().get(i).getTraceID();
            double currentAImpact = a.getTraceVector().get(i).getImpact();

            a_traceVector.add(new TraceTuple(currentAID, influence_factor_old * currentAImpact));
            i++;
        }

        a = new TraceableBoolean(new_a_value, a_traceVector);
    }

    /** Destructively mutates the individual in some default manner.  The default form
     does a bit-flip with a probability depending on parameters. */
    public void defaultMutate(EvolutionState state, int thread)
    {
        BitVectorSpecies s = (BitVectorSpecies)species;  // where my default info is stored
        for(int x=0;x<genome.length;x++)
        {
            if (state.random[thread].nextBoolean(s.mutationProbability(x)))
            {
                TraceableBoolean old = genome[x];
                for(int retries = 0; retries < s.duplicateRetries(x) + 1; retries++)
                {
                    switch(s.mutationType(x))
                    {
                        case BitVectorSpecies.C_FLIP_MUTATION:
                            state.mutationCounter--;
                            mutateGene(genome[x], !genome[x].getValue(), state.mutationCounter);
                            break;
                        case BitVectorSpecies.C_RESET_MUTATION:
                            state.mutationCounter--;
                            genome[x] = new TraceableBoolean(state.random[thread].nextBoolean(), state.mutationCounter);
                            break;
                        default:
                            state.output.fatal("In TracableBitVectorIndividual.defaultMutate, default case occurred when it shouldn't have");
                            break;
                    }
                    if (genome[x].getValue() != old.getValue()) break; //NOTE: i think the values are needed here
                    // else genome[x] = old;  // try again
                }
            }
        }
    }

    /** Initializes the individual by randomly flipping the bits */
    public void reset(EvolutionState state, int thread)
    {
        for(int x=0;x<genome.length;x++) {
            genome[x] = new TraceableBoolean(state.random[thread].nextBoolean(), this.traceID);
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
        s.append( Code.encode( genome.length ) );
        for( int i = 0 ; i < genome.length ; i++ ) {
            s.append(Code.encode(genome[i].toString()));
        }
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
        for( int i = 0 ; i < genome.length ; i++ )
            s.append( Code.encode( genome[i].toString()));
        return s.toString();
    }

    protected void parseGenotype(final EvolutionState state,
                                 final LineNumberReader reader) throws IOException
    {
        String s = reader.readLine();
        DecodeReturn d = new DecodeReturn(s);
        Code.decode( d );

        int lll = (int)(d.l);

        genome = new TraceableBoolean[ lll ];

        // read in the genes
        for( int i = 0 ; i < genome.length ; i++ )
        {
            Code.decode( d );
            genome[i] = new TraceableBoolean();
            genome[i].fromString(d.s);
        }
    }

    public boolean equals(Object ind)
    {
        System.out.println("equals called"); //DEBUGG: check if called, delete if i know when
        if (ind==null) return false;
        if (!(this.getClass().equals(ind.getClass()))) return false; // SimpleRuleIndividuals are special.
        TraceableBitVectorIndividual i = (TraceableBitVectorIndividual)ind;
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
    { genome = (TraceableBoolean[]) gen; }
    public int genomeLength()
    { return genome.length; }

    public void setGenomeLength(int len)
    {
        System.out.println("setGenotypeLength called"); //DEBUGG: check if called, delete if i know when
        TraceableBoolean[] newGenome = new TraceableBoolean[len];
        System.arraycopy(genome, 0, newGenome, 0,
                genome.length < newGenome.length ? genome.length : newGenome.length);
        genome = newGenome;
    }

    public void writeGenotype(final EvolutionState state,
                              final DataOutput dataOutput) throws IOException
    {
        System.out.println("writeGenotype called"); //DEBUGG: check if called, delete if i know when
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
        genome = new TracableBoolean[len];
    for(int x=0;x<genome.length;x++)
        genome[x] = dataInput.readBoolean();
     */
    }

    /** Implements distance as hamming distance. */
    public double distanceTo(Individual otherInd)
    {
        if (!(otherInd instanceof TraceableBitVectorIndividual))
            return super.distanceTo(otherInd);  // will return infinity!

        TraceableBitVectorIndividual other = (TraceableBitVectorIndividual) otherInd;
        TraceableBoolean[] otherGenome = other.genome;
        double hammingDistance =0;
        for(int i=0; i < other.genomeLength(); i++)
        {
            if(genome[i].getValue() ^ otherGenome[i].getValue())  //^ is xor //TODO: was this right?
                hammingDistance++;
        }

        return hammingDistance;
    }
}
