/*
  Copyright 2006 by Sean Luke
  Licensed under the Academic Free License version 3.0
  See the file "LICENSE" for more information
*/


package ec.eda;
import ec.*;
import ec.simple.*;

/* 
 * CMAESInitializer.java
 * 
 * Created: Wed Jul  8 12:35:31 EDT 2015
 * By: Sam McKay and Sean Luke
 */

/**
 * SimpleInitializer is a default Initializer which initializes a Population
 * by calling the Population's populate(...) method.  For most applications,
 * this should suffice.
 *
 * @author Sean Luke
 * @version 1.0 
 */

public class CMAESInitializer extends SimpleInitializer
    {
    private static final long serialVersionUID = 1;

    public Population setupPopulation(final EvolutionState state, int thread)
        {
        Population p = super.setupPopulation(state, thread);
        
        // reset to lambda in size!
        for(int i = 0; i < p.subpops.length; i++)
            if (p.subpops[i].species instanceof CMAESSpecies)
                p.subpops[i].individuals = new Individual[(int)(((CMAESSpecies)p.subpops[i].species).lambda)];
                                        
        return p;
        }
    }
