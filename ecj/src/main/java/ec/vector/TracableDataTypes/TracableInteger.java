package ec.vector.TracableDataTypes;

public class TracableInteger{

    public TracableInteger(){
        _value = 0;
        _traceID = 0;
    }

    public TracableInteger(int value, int traceID){
        _value = value;
        _traceID = traceID;
    }

    private int _value;
    private int _traceID;

    public int getValue(){ return _value; }
    public int getTraceID(){ return _traceID; }

    public void setValue(int value, int traceID){ this._value = value; this._traceID = traceID; }

}
