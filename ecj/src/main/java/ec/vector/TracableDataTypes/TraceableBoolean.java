package ec.vector.TracableDataTypes;

import java.util.ArrayList;
import java.util.List;

public class TraceableBoolean {

    public TraceableBoolean(){
        _value = false;
        _traceVector = new ArrayList<TraceTuple>();
    }

    public TraceableBoolean(boolean value, int initialTraceID){
        _value = value;
        _traceVector = new ArrayList<TraceTuple>();
        _traceVector.add(new TraceTuple(initialTraceID, 1.0));
    }

    public TraceableBoolean(boolean value, List<TraceTuple> traceIDs){
        _value = value;
        _traceVector = new ArrayList<TraceTuple>();
        for(TraceTuple tt : traceIDs){
            _traceVector.add(tt);
        }
    }

    private boolean _value;
    private List<TraceTuple> _traceVector;

    public boolean getValue(){ return _value; }
    public List<TraceTuple> getTraceVector(){ return _traceVector; }

    public void setValue(boolean value, List<TraceTuple> traceVector){ this._value = value; this._traceVector = traceVector; }

}
