package ec.vector.TracableDataTypes;

import java.util.ArrayList;
import java.util.List;



public class TraceableFloat{

    public TraceableFloat(){
        _value = 0;
        _traceVector = new ArrayList<TraceTuple>();
    }

    public TraceableFloat(float value, int initialTraceID){
        _value = value;
        _traceVector = new ArrayList<TraceTuple>();
        _traceVector.add(new TraceTuple(initialTraceID, 1.0));
    }

    public TraceableFloat(float value, List<TraceTuple> traceIDs){
        _value = value;
        _traceVector = new ArrayList<TraceTuple>();
        for(TraceTuple tt : traceIDs){
            _traceVector.add(tt);
        }
    }

    private float _value;
    private List<TraceTuple> _traceVector;

    public float getValue(){ return _value; }
    public List<TraceTuple> getTraceVector(){ return _traceVector; }

    public void setValue(float value, List<TraceTuple> traceVector){ this._value = value; this._traceVector = traceVector; }

}
