package ec.vector.TracableDataTypes;

import java.util.ArrayList;
import java.util.List;



public class TraceableString {

    public TraceableString(){
        _value = "";
        _traceVector = new ArrayList<TraceTuple>();
    }

    public TraceableString(String value, int initialTraceID){
        _value = value;
        _traceVector = new ArrayList<TraceTuple>();
        _traceVector.add(new TraceTuple(initialTraceID, 1.0));
    }

    public TraceableString(String value, List<TraceTuple> traceIDs){
        _value = value;
        _traceVector = new ArrayList<TraceTuple>();
        for(TraceTuple tt : traceIDs){
            _traceVector.add(tt);
        }
    }

    private String _value;
    private ArrayList<TraceTuple> _traceVector;

    public String getValue(){ return _value; }
    public ArrayList<TraceTuple> getTraceVector(){ return _traceVector; }

    public void setValue(String value, ArrayList<TraceTuple> traceVector){ this._value = value; this._traceVector = traceVector; }



    public String toString(){
        String outString = "";
        outString += "["+_value+",";
        for(TraceTuple tt : _traceVector){
            outString += tt.toString();
        }
        outString += "]";
        return outString;
    }

    public void fromString(String data){

        //TODO: check if the input string has the right format!

        String splitter = ",\\["; //otherwise you get a unclosed brackets error...

        //remove the first and last character and split the string
        String[] parsedData = data.substring(1, data.length() - 2).split(splitter);
        _value = parsedData[0];
        _traceVector = new ArrayList<TraceTuple>();
        for(int i = 1; i < parsedData.length; i++){
            TraceTuple toAdd = new TraceTuple(0, 1.0); //these are dummy values which are overwritten in the next line
            toAdd.fromString("["+parsedData[i]); //re add the "["
            _traceVector.add(toAdd);
        }
    }
}
