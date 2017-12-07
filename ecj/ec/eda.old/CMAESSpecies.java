/*
  Copyright 2015 by Sean Luke
  Licensed under the Academic Free License version 3.0
  See the file "LICENSE" for more information
*/


package ec.eda;
import ec.*;
import ec.vector.*;
import ec.util.*;
import ec.simple.SimpleFitness;
import java.io.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.ops.CommonOps;
import org.ejml.ops.NormOps;
import org.ejml.ops.RandomMatrices;
import org.ejml.simple.SimpleMatrix;

import java.util.Arrays;

/* 
 * CMAESSpecies.java
 * 
 * Created: Wed Jul  8 12:29:50 EDT 2015
 * By: Sam McKay and Sean Luke
 */

/**
 * Species is a prototype which defines the features for a set of individuals
 * in the population.  Typically, individuals may breed if they belong to the
 * same species (but it's not a hard-and-fast rule).  Each Subpopulation has
 * one Species object which defines the species for individuals in that
 * Subpopulation.
 *
 * <p>Species are generally responsible for creating individuals, through
 * their newIndividual(...) method.  This method usually clones its prototypical
 * individual and makes some additional modifications to the clone, then returns it.
 * Note that the prototypical individual does <b>not need to be a complete individual</b> --
 * for example, GPSpecies holds a GPIndividual which doesn't have any trees (the tree
 * roots are null).
 *
 * <p>Species also holds a prototypical breeding pipeline meant to breed
 * this individual.  To breed individuals of this species, clone the pipeline
 * and use the clone.

 <p><b>Parameters</b><br>
 <table>
 <tr><td valign=top><i>base</i>.<tt>ind</tt><br>
 <font size=-1>classname, inherits and != ec.Individual</font></td>
 <td valign=top>(the class for the prototypical individual for the species)</td></tr>

 <tr><td valign=top><i>base</i>.<tt>fitness</tt><br>
 <font size=-1>classname, inherits and != ec.Fitness</font></td>
 <td valign=top>(the class for the prototypical fitness for the species)</td></tr>

 <tr><td valign=top><i>base</i>.<tt>numpipes</tt><br>
 <font size=-1>int &gt;= 1</font></td>
 <td valign=top>(total number of breeding pipelines for the species)</td></tr>

 <tr><td valign=top><i>base</i>.<tt>pipe</tt><br>
 <font size=-1>classname, inherits and != ec.BreedingPipeline</font></td>
 <td valign=top>(the class for the prototypical Breeding Pipeline)</td></tr>

 </table>


 <p><b>Parameter bases</b><br>
 <table>
 <tr><td valign=top><i>base</i>.<tt>ind</tt></td>
 <td>i_prototype (the prototypical individual)</td></tr>

 <tr><td valign=top><i>base</i>.<tt>pipe</tt></td>
 <td>pipe_prototype (breeding pipeline prototype)</td></tr>

 <tr><td valign=top><i>base</i>.<tt>fitness</tt></td>
 <td>f_prototype (the prototypical fitness)</td></tr>

 </table>



 * @author Sam McKay and Sean Luke
 * @version 1.0 
 */

public class CMAESSpecies extends FloatVectorSpecies
    {
    public static final String P_SPECIES = "species";
    
    public static final String P_LAMBDA = "lambda";
    public static final String P_MU = "mu";
    public static final String P_SIGMA = "sigma";
    public static final String P_MEAN = "mean";
    public static final String P_WEIGHTS = "weight";

    public static final String P_CC = "cc";
    public static final String P_CS = "cs";
    public static final String P_C1 = "c1";
    public static final String P_CMU = "cmu";
    public static final String P_DAMPS = "damps";
    
    // SEAN: All of this desperately asks for documentation and explanation

    public SimpleMatrix b;
    public SimpleMatrix d;
    public SimpleMatrix c;
    public DenseMatrix64F bd;
    public DenseMatrix64F sbd;
    public SimpleMatrix invsqrtC;
    public SimpleMatrix invD;

    public SimpleMatrix xmean;
    public double sigma;

    public double lambda;
    public int mu;
    public double[] weights;
    public double mueff;

    // Strategy parameter setting: Adaptation
    public double cc;  // time constant for cumulation for C
    public double cs;  // t-const for cumulation for sigma control
    public double c1;  // learning rate for rank-one update of C
    public double cmu;  // and for rank-mu update
    public double damps; // damping for sigma

    public SimpleMatrix ps;
    public SimpleMatrix pc;
    public double chiN;

    // SEAN: you can replace counteval with state.generation
    // SAM: more correctly, (int)(state.generation * lambda) if it were to match the previous code
    // However, things can be simplified to just using the generation number
    // I'm not sure why it didn't do that in the first place

    // public int counteval;
    public int eigeneval;

    public Parameter defaultBase()
        {
        return EDADefaults.base().push(P_SPECIES);
        }

    public void setup(final EvolutionState state, final Parameter base)
        {
        super.setup(state, base);
        MersenneTwisterFast random = state.random[0];

        Parameter def = defaultBase();

        Parameter subpopBase = base.pop();
        Parameter subpopDefaultBase =  ECDefaults.base().push(Subpopulation.P_SUBPOPULATION);

        // set myself up and define my initial distribution here
        int n = genomeSize;
        b = SimpleMatrix.identity(n);
        c = SimpleMatrix.identity(n);
        d = SimpleMatrix.identity(n);
        bd = CommonOps.identity(n,n);
        sbd = CommonOps.identity(n,n);
        invsqrtC = SimpleMatrix.identity(n);
        invD = SimpleMatrix.identity(n);

        // Initialize dynamic (internal) strategy parameters and constants
        pc = new SimpleMatrix(n,1);
        ps = new SimpleMatrix(n,1);   // evolution paths for C and sigma
        chiN=Math.pow(n,0.5)*(1.0-1.0/(4.0*n)+1.0/(21.0*n*n));  // expectation of ||N(0,I)|| == norm(randn(N,1))

        xmean = new SimpleMatrix(genomeSize,1);
        
        // is there a better way to initialize the mean?
        if(state.parameters.exists(base.push(P_MEAN).push(""+0),def.push(P_MEAN).push(""+0)))
            {
            try 
                {
                for(int i = 0; i < genomeSize; i++)
                    {
                    double m_i = state.parameters.getDouble(base.push(P_MEAN).push(""+i),def.push(P_MEAN).push(""+i));
                    xmean.set(i,0,m_i);
                    }
                // SEAN: what if this parameter isn't a number?  You need to throw an error like I discussed in the email
                // SAM: Changed, although I'm not sure this is the correct to to handle it.
                }
            catch (NumberFormatException e) 
                {
                state.output.fatal("Invalid value for CMA-ES mean: \n" + e);
                }
            } 
        else
            {
            // SEAN: Why are you using random.nextDouble()?  This is a strange average
            // location for initial means (<0.5, 0.5, ... >).  I'd just make it all <0,0,0,0> by default.
            // Also I'd issue a warning (state.output.warning(...)) that the mean is missing and so it'll
            // be set to 0.

            // SAM: that was to correspond to Matlab code. Changed

            state.output.warning("CMA-ES mean not provided, using all zeros", base.push(P_MEAN).push("" + 0), def.push(P_MEAN).push("" + 0));
            // xmean was 0 initialized previously
            }

        // SEAN: what if this parameter is missing?  You're going to get an error thrown
        // here.  Instead say [ugh]
        // if (!state.parameters.exists(base.push(P_SIGMA),def.push(P_SIGMA))
        //     issue a warning, say sigma will be set to (I dunno) 1.0
        // else
        //     sigma = state.parameters.getDouble(base.push(P_SIGMA),def.push(P_SIGMA), 0.0);
        //     if (sigma == 0) state that sigma can't be 0.0, issue a fatal
        //     if (sigma < 0) state that sigma is not a number, issue a fatal 
        
        // SAM: changed

        if(!state.parameters.exists(base.push(P_SIGMA),def.push(P_SIGMA)))
            {
            state.output.warning("CMA-ES sigma was not provided, defaulting to value of 1.0",base.push(P_SIGMA),def.push(P_SIGMA));
            sigma = 1.0;
            }
        else
            {
            sigma = state.parameters.getDouble(base.push(P_SIGMA),def.push(P_SIGMA),0.0);
            if (sigma == 0) 
                state.output.fatal("CMA-ES sigma may not be set to 0",base.push(P_SIGMA),def.push(P_SIGMA));
            if (sigma < 0) 
                state.output.fatal("CMA-ES sigma is not a valid number",base.push(P_SIGMA),def.push(P_SIGMA));
            }


        // SEAN: Subpopulations haven't been set up yet.  Since CMA-ES has a default
        // setting for lambda, we should use that unless it's been specified by the user.
        // So this brings up an interesting question: what do we do about Subpopulation.P_SUBPOPSIZE?
        // Normally an evolution strategy in ECJ uses P_SUBPOPSIZE to set the *initial* subpopulation
        // size, and uses "lambda" to set the *later* subpopulation sizes.  This isn't relevant
        // for CMA-ES.  We could do it in a few ways:
        // 1. Ignore P_SUBPOPSIZE and resize the subpopulation individuals array as necessary in 
        //    CMAESBreeder.breedPopulation(...).  Issue a warning reminding the user that this is the
        //    case.  The danger here is that if someone relies on the previous subpopulation size, they'll get broken.
        // 2. Set the P_SUBPOPSIZE parameter (using parameters.set(...)) to our desired lambda,
        //    issuing a warning of course.  This feels hokey.
        // 3. Use P_SUBPOPSIZE for our lambda.  This prevents us from using the CMA-ES default.
        //
        // I would do  #1.  At present the only two classes which rely on P_SUBPOPSIZE other
        // than Subpopulation are MuCommaLambdaBreeder (which looks it up to set lambda to it
        // if the user didn't specify lambda), and NSGAEvaluator.  So I think we're safe.
        // 
        // So do this: if lambda isn't specified as P_LAMBDA, then use the CMA-ES default.
        // otherwise use what the user specified.

        // SAM: changed here and in CMAESBreeder.java to use #1
        // THIS DIDN'T WORK. When initially populated, the subpopulation only has he amount specified by P_SUBPOPSIZE

        if(!state.parameters.exists(base.push(P_LAMBDA),def.push(P_LAMBDA)))
            {
            // SAM: should there be a warning here?
            state.output.warning("CMA-ES lambda was not provided, resorting to default to value",base.push(P_LAMBDA),def.push(P_LAMBDA));
            lambda = 4+Math.floor(3*Math.log(n));
            }
        else
            {
            lambda = state.parameters.getDouble(base.push(P_LAMBDA),def.push(P_LAMBDA),0.0);
            if (lambda <= 0) 
                state.output.fatal("If the CMA-ES lambda parameter is provided, it must be a valid number, and not 0.0",base.push(P_LAMBDA),def.push(P_LAMBDA));
            }
            
        // SEAN: Don't do this, what if the parameter is broken or not a number?  I'd
        // test to see if it's missing, and use lambda/2.0 if so.  BTW, it should be Math.floor(lambda/2.0), right?

        // SAM: changed

        if(!state.parameters.exists(base.push(P_MU),def.push(P_MU)))
            {
            // SAM: should there be a warning here?
            state.output.warning("CMA-ES mu was not provided, resorting to default to value",base.push(P_MU),def.push(P_MU));
            mu = (int)(Math.floor(lambda/2.0));
            }
        else
            {
            mu = state.parameters.getInt(base.push(P_MU),def.push(P_MU),0);
            if (mu <= 0) 
                state.output.fatal("If the CMA-ES mu parameter is provided, it must be a valid integer, and not 0",base.push(P_MU),def.push(P_MU));
            }

        // SEAN: why use Math.floor(mu) here and elsewhere.  Isn't mu supposed to be already floored?
                
        // SAM: That was the way it worked in the Matlab code. When I was check my code's values against it 
        // the Math.log(mu+1.0/2.0) produced slightly different values.
        
        // SAM: Changed, mu is now an integer

        // should we be able to set custom weights?
        weights = new double[mu];
        double sum = 0.0;
        double sumSqr = 0.0;
        if(state.parameters.exists(base.push(P_WEIGHTS).push(""+0),def.push(P_WEIGHTS).push(""+0)))
            {
            try 
                {
                for(int i = 0; i < genomeSize; i++)
                    {
                    weights[i] = state.parameters.getDouble(base.push(P_WEIGHTS).push(""+i),def.push(P_WEIGHTS).push(""+i));
                    sum += weights[i];
                    }
                }
            catch (NumberFormatException e) 
                {
                state.output.fatal("Invalid value for CMA-ES weight: \n" + e);
                }
            }
        else
            {
            for(int i = 0; i < mu; i++)
                {
                weights[i] = Math.log(mu+1.0/2.0)-Math.log(i+1.0);
                sum += weights[i];
                }
            for(int i = 0; i < mu; i++)
                {
                weights[i] /= sum;
                sumSqr += weights[i]*weights[i];
                }
            }



        // SEAN: ????
        // SAM: removed the mu = Math.floor(mu) that I assume was the reason for te question marks
        mueff=1.0/sumSqr;

        // SAM: changed to allow initialization of the following parameters
        // should these values be checked in any way to make sure they make sense?

        // Strategy parameter setting: Adaptation

        if(!state.parameters.exists(base.push(P_CC),def.push(P_CC)))
            {
            cc = (4.0+mueff/n) / (n+4.0 + 2.0*mueff/n);  // time constant for cumulation for C
            }
        else
            {
            cc = state.parameters.getDouble(base.push(P_CC),def.push(P_CC),0.0);
            if (cc <= 0) 
                state.output.fatal("If the CMA-ES cc parameter is provided, it must be a valid number, and not 0.0",base.push(P_CC),def.push(P_CC));
            }

        if(!state.parameters.exists(base.push(P_CS),def.push(P_CS)))
            {
            cs = (mueff+2.0)/(n+mueff+5.0);  // t-const for cumulation for sigma control
            }
        else
            {
            cs = state.parameters.getDouble(base.push(P_CS),def.push(P_CS),0.0);
            if (cs <= 0) 
                state.output.fatal("If the CMA-ES cs parameter is provided, it must be a valid number, and not 0.0",base.push(P_CS),def.push(P_CS));
            }

        if(!state.parameters.exists(base.push(P_C1),def.push(P_C1)))
            {
            c1 = 2.0 / ((n+1.3)*(n+1.3)+mueff);  // learning rate for rank-one update of C
            }
        else
            {
            c1 = state.parameters.getDouble(base.push(P_C1),def.push(P_C1),0.0);
            if (c1 <= 0) 
                state.output.fatal("If the CMA-ES c1 parameter is provided, it must be a valid number, and not 0.0",base.push(P_C1),def.push(P_C1));
            }
        
        if(!state.parameters.exists(base.push(P_CMU),def.push(P_CMU)))
            {
            cmu = Math.min(1.0-c1, 2.0*(mueff-2.0+1.0/mueff) / ((n+2.0)*(n+2.0)+mueff));
            }
        else
            {
            cmu = state.parameters.getDouble(base.push(P_CMU),def.push(P_CMU),0.0);
            if (cmu <= 0) 
                state.output.fatal("If the CMA-ES cmu parameter is provided, it must be a valid number, and not 0.0",base.push(P_CMU),def.push(P_CMU));
            }

        if(!state.parameters.exists(base.push(P_DAMPS),def.push(P_DAMPS)))
            {
            damps = 1.0 + 2.0*Math.max(0.0, Math.sqrt((mueff-1.0)/(n+1.0))-1.0) + cs; // damping for sigma
            }
        else
            {
            damps = state.parameters.getDouble(base.push(P_DAMPS),def.push(P_DAMPS),0.0);
            if (damps <= 0) 
                state.output.fatal("If the CMA-ES damps parameter is provided, it must be a valid number, and not 0.0",base.push(P_DAMPS),def.push(P_DAMPS));
            }

        }

    public Object clone()
        {
        CMAESSpecies myobj = (CMAESSpecies) (super.clone());
            
        // clone the distribution and other variables here
        myobj.c = c.copy();
        myobj.b = b.copy();
        myobj.d = d.copy();
        myobj.bd = bd.copy();
        myobj.sbd = sbd.copy();
        myobj.invsqrtC = invsqrtC.copy();
        myobj.invD = invD.copy();
    
        myobj.xmean = xmean.copy();
        myobj.ps = ps.copy();
        myobj.pc = pc.copy();

        // SEAN: there is no need to say this stuff.  Clone already copies over pointer
        // references.  Also, did you mean to not copy xmean?
        // SAM: Changed, and that was a mistake, although I don't believe it would matter
    

        // myobj.sigma = sigma;

        // myobj.lambda = lambda;
        // myobj.mu = mu;
        
        // SEAN: do you expect weights to be changed?  I think they shouldn't be cloned here.
        // SAM: I haven't been making that assumption, but that is almost certainly the case.


        // myobj.weights = weights.clone();
        // myobj.mueff = mueff;

        // myobj.cc = cc;
        // myobj.cs = cs;
        // myobj.c1 = c1;
        // myobj.cmu = cmu;
        // myobj.damps = damps;

        // myobj.chiN = chiN;

        // myobj.counteval = counteval;
        // myobj.eigeneval = eigeneval;
            
        return myobj;
        } 

    public Individual newIndividual(final EvolutionState state, int thread)
        {
        Individual newind = super.newIndividual(state, thread);
        MersenneTwisterFast random = state.random[thread];

        // need error check? Should be checked at the beginning
        // if (!(newind instanceof DoubleVectorIndividual))  // uh oh
        //     state.output.fatal("To use CMAESSpecies, the species must be initialized with a DoubleVectorIndividual.  But it contains a " + newind);
        
        DoubleVectorIndividual dvind = (DoubleVectorIndividual)(newind);

        DenseMatrix64F genome = DenseMatrix64F.wrap(genomeSize,1,dvind.genome);
        DenseMatrix64F temp = new DenseMatrix64F(genomeSize,1);

        // arz(:,k) = randn(N,1); % standard normally distributed vector
        // arx(:,k) = xmean + sigma*(B*D*arz(:,k));
        while(true)
            {           

            for( int i = 0; i < genomeSize; i++ ) 
                dvind.genome[i] = random.nextGaussian();

            CommonOps.mult(sbd,genome,temp); // temp = sigma*b*d*genome;
            CommonOps.add(temp,xmean.getMatrix(),genome); // genome = temp + xmean;

            for(int i = 0; i < genomeSize; i++)
                if(dvind.genome[i] < minGene[i] || dvind.genome[i] > maxGene[i]) continue;

            return newind;
            }
        }

    public void updateDistribution(final EvolutionState state, final Subpopulation subpop)
        {
        // % Sort by fitness and compute weighted mean into xmean
        // [arfitness, arindex] = sort(arfitness); % minimization
        // xmean = arx(:,arindex(1:mu))*weights;   % recombination            % Eq.39
        // counteval += lambda;

        // only need partial sort?
        Arrays.sort(subpop.individuals);

        SimpleMatrix artmp = new SimpleMatrix(genomeSize,(int)(mu));
        SimpleMatrix xold = xmean;
        xmean = new SimpleMatrix(genomeSize, 1);

        for(int i = 0; i < mu; i++)
            {
            DoubleVectorIndividual dvind = (DoubleVectorIndividual)(subpop.individuals[i]);
            SimpleMatrix arz = new SimpleMatrix(DenseMatrix64F.wrap(genomeSize,1,dvind.genome));
            // arz.transpose().print("%25.15e");
            arz = (arz.minus(xold).divide(sigma));
            for(int j = 0; j < genomeSize; j++) 
                {
                xmean.set(j,0, xmean.get(j,0)+weights[i]*dvind.genome[j]);
                artmp.set(j,i,arz.get(j,0));
                }
            }
            
        // SEAN: In general this math is pretty hard to read.
        // SAM: I can break it up, but I don't think it's going to be easier to read

        // % Cumulation: Update evolution paths
        ps = ps.scale(1.0-cs).plus( invsqrtC.mult( xmean.minus(xold).scale(1.0/sigma) ).scale( Math.sqrt( cs*(2.0-cs) * mueff ) ) );
        int hsig = (ps.elementMult(ps).elementSum() / (1.0-Math.pow(1.0-cs,2.0*state.generation))/genomeSize < 2.0+4.0/(genomeSize+1))?1:0;
        pc = pc.scale(1.0-cc).plus((xmean.minus(xold).scale(1.0/sigma)).scale(hsig).scale(Math.sqrt(cc*(2.0-cc)*mueff)));

        // % Adapt covariance matrix C
        c = c.scale(1.0-c1-cmu).plus(
            pc.mult(pc.transpose()).plus(c.scale(1.0-hsig).scale(cc).scale(2.0-cc)).scale(c1)).plus(
                (artmp.mult(SimpleMatrix.diag(weights).mult(artmp.transpose()))).scale(cmu));

        // % Adapt step-size sigma
        sigma = sigma*Math.exp((cs/damps)*(ps.normF()/chiN - 1.0));


        // % Update B and D from C
        // SEAN: don't do x/y/z/d, it's difficult to identify order of operations in division
        // SAM: Changed to x/(y*z*d)
        if((state.generation - eigeneval) > 1.0/((c1+cmu)*genomeSize*10.0) )
            {
            eigeneval = state.generation;

            // make sure the matrix is symmetric (it should be already)
            // not sure if this is necessary           
            for(int i = 0; i < genomeSize; i++)
                for(int j = 0; j < i; j++)
                    c.set(j,i,c.get(i,j));

            // this copy gets modified by the decomposition
            DenseMatrix64F copy = c.copy().getMatrix();
            EigenDecomposition<DenseMatrix64F> eig = DecompositionFactory.eig(genomeSize,true,true);
            eig.decompose(copy);

            SimpleMatrix dinv = new SimpleMatrix(genomeSize,genomeSize);
            for(int i = 0; i < genomeSize; i++)
                {
                double eigrt = Math.sqrt(eig.getEigenvalue(i).real);
                d.set(i,i,eigrt);
                dinv.set(i,i,1/eigrt);
                CommonOps.insert(eig.getEigenVector(i),b.getMatrix(),0,i);
                }

            invsqrtC = b.mult(dinv.mult(b.transpose()));
            CommonOps.mult(b.getMatrix(),d.getMatrix(),bd);
            }

        CommonOps.scale(sigma, bd, sbd);

        // % Break, if fitness is good enough or condition exceeds 1e14, better termination methods are advisable 
        // if arfitness(1) <= stopfitness || max(D) > 1e7 * min(D)
        //   break;
        // end
        if(CommonOps.elementMax(d.extractDiag().getMatrix()) > 1e7*CommonOps.elementMin(d.extractDiag().getMatrix()))
            {
            state.evaluator.setRunComplete("CMAESSpecies: Stopped because matrix condition exceeded limit.");
            }
        }



    }

