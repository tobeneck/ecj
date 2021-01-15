package ec.vector.TracableDataTypes;

import java.util.ArrayList;
import java.util.List;

public class TraceableInteger {

    public TraceableInteger(){
        _value = 0;
        _traceVector = new ArrayList<TraceTuple>();
    }

    public TraceableInteger(int value, int initialTraceID){
        _value = value;
        _traceVector = new ArrayList<TraceTuple>();
        _traceVector.add(new TraceTuple(initialTraceID, 1.0));
    }

    public TraceableInteger(int value, List<TraceTuple> traceIDs){
        _value = value;
        _traceVector = new ArrayList<TraceTuple>();
        for(TraceTuple tt : traceIDs){
            _traceVector.add(tt);
        }
    }

    private int _value;
    private List<TraceTuple> _traceVector;

    public int getValue(){ return _value; }
    public List<TraceTuple> getTraceVector(){ return _traceVector; }

    public void setValue(int value, List<TraceTuple> traceID){ this._value = value; this._traceVector = traceID; }

}
