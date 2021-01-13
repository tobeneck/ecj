package ec.vector.TracableDataTypes;

public class TraceTuple{
    public TraceTuple(int traceID, double impact){
        _traceID = traceID;
        _impact = impact;
    }

    private int _traceID;
    private double _impact;

    public void setTraceID(int traceID){ _traceID = traceID; }
    public void setImpact(double impact){ _impact = impact; }

    public int getTraceID(){ return _traceID; }
    public double getImpact(){ return _impact; }
}
