package ec.vector.TracableDataTypes;

public class TracableBoolean{

    public TracableBoolean(){
        _value = false;
        _traceID = 0;
    }

    public TracableBoolean(boolean value, int traceID){
        _value = value;
        _traceID = traceID;
    }

    private boolean _value;
    private int _traceID;

    public boolean getValue(){ return _value; }
    public int getTraceID(){ return _traceID; }

    public void setValue(boolean value, int traceID){ this._value = value; this._traceID = traceID; }

}
