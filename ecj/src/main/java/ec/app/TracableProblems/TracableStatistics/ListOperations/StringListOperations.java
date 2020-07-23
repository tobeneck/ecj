package ec.app.TracableProblems.TracableStatistics.ListOperations;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class StringListOperations {

    /**
     * own implementation
     * Returns the entropy Shannon Entropy of a given StringList
     */
    public static Double getShannonEntropy(List<String> strList){
        List<Integer> intList = new ArrayList<Integer>();
        for(String s : strList)
            intList.add(Integer.valueOf(s));

        Collections.sort(intList);

        int currentK = intList.get(0);
        double currentP = 1.0/intList.size();
        double sum = 0.0;

        for(int i = 1; i < intList.size(); i++){
            if(currentK == intList.get(i))
                currentP += 1.0/intList.size();
            else {
                sum += (-1.0) * currentP * DoubleListOperations.log2(currentP);
                currentP = 1.0/intList.size(); //since it is now already observed once
                currentK = intList.get(i);
            }
        }

        return sum;
    }
}
