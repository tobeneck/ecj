package ec.app.TracableProblems.TracableStatistics.ListOperations;

import java.util.ArrayList;
import java.util.List;

public class DoubleListOperations {
    /**
     * returns the log2 if a double
     */
    public static double log2(double a) {
        return Math.log(a) / Math.log(2);
    }


    /**
     * normalizes a list so the overall value is 1
     */
    public static List<Double> normalizeList(List<Double> inputList){
        if(getMin(inputList) < 0)
            throw new IndexOutOfBoundsException("Impact List has a element < 0! " + getMin(inputList));


        double sum = getSum(inputList);

        List<Double> returnList = new ArrayList<Double>();
        for(double d : inputList)
            returnList.add(d / sum);

        return returnList;
    }


    public static double getSum(List<Double> inputList){
        double sum = 0.0;
        for(double d : inputList)
            sum += d;
        return sum;
    }

    /**
     * returns the minimum value of a double List
     */
    public static double getMin(List<Double> inputList){
        double min = inputList.get(0);

        for(int i = 1; i < inputList.size(); i++)
            min = Math.min(inputList.get(i),min);

        return min;
    }

    /**
     * returns the maximum value of a double List
     */
    public static double getMax(List<Double> inputList){
        double max = inputList.get(0);

        for(int i = 1; i < inputList.size(); i++)
            max = Math.max(inputList.get(i),max);

        return max;
    }

    /**
     * returns the List as a CSV String
     */
    public static String listToString(List<Double> inputList){
        String returnString = "";

        for(double d : inputList)
            if(returnString.isEmpty())
                returnString += d;
            else
                returnString += ", " + d;

        return returnString;
    }



}
