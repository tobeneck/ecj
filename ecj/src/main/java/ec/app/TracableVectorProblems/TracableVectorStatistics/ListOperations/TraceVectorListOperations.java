package ec.app.TracableVectorProblems.TracableVectorStatistics.ListOperations;

import ec.vector.TracableDataTypes.TraceTuple;

import java.util.ArrayList;
import java.util.List;

public class TraceVectorListOperations {

    /**
     * Compares if two traceVectors are equal
     * @param traceVector1 The first vector to be compared
     * @param traceVector2 The second vector to be compared
     * @return vector1 and vector2 equal?
     */
    private static boolean areTraceVectorsEqual(List<TraceTuple> traceVector1, List<TraceTuple> traceVector2){
        if(traceVector1.size() != traceVector2.size())
            return false;

        for(int i = 0; i < traceVector1.size(); i++){
            if(traceVector1.get(i).getTraceID() != traceVector2.get(i).getTraceID()
                    || traceVector1.get(i).getImpact() != traceVector2.get(i).getImpact())
                return false;
        }

        return true;
    }

    /**
     * returns of a specific
     * @param traceVectors the traceVectors to be checked
     * @param traceVector the traceVector to search for
     * @return
     */
    private static boolean traceVectorsContainTraceVector(ArrayList<ArrayList<TraceTuple>> traceVectors, ArrayList<TraceTuple> traceVector){
        for(int i = 0; i < traceVectors.size(); i++){
            if(areTraceVectorsEqual(traceVectors.get(i), traceVector))
                return true;
        }
        return false;
    }

    /**
     * own implementation
     * Returns the entropy Shannon Entropy of a set of traceVectors
     * @param traceVectors the traceVectors of one genome. NOT the traceVectors of one individual!
     * @return shannon entropy of the genome
     */
    public static Double getShannonEntropy(ArrayList<ArrayList<TraceTuple>> traceVectors){

        //compute the alphabet, e.g. a vector only containing every different traceVector once:
        ArrayList<ArrayList<TraceTuple>> alphabet = new ArrayList<ArrayList<TraceTuple>>();
        for(int i = 0; i < traceVectors.size(); i++){
            ArrayList<TraceTuple> currentTraceVector = traceVectors.get(i);
            if(!traceVectorsContainTraceVector(alphabet, currentTraceVector))
                alphabet.add(currentTraceVector);
        }

        double currentP = 0.0;//1.0/intList.size();
        double sum = 0.0;


        for(int k = 0; k < alphabet.size(); k++){//iterate over the alphabet to compute for the current gene
            List<TraceTuple> currentLetter = alphabet.get(k);
            for(int i = 0; i < traceVectors.size(); i++){//again iterate over the genome to compute the current P
                if(areTraceVectorsEqual(currentLetter, traceVectors.get(i))) //
                    currentP += 1.0/(double)traceVectors.size();
            }

            sum += (-1.0) * currentP * DoubleListOperations.log2(currentP);
            currentP = 0.0;
        }
        return sum;
    }
}
