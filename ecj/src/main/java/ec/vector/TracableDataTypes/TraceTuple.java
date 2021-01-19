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

    public String toString(){
        return "["+_traceID+","+_impact+"]";
    }


    /**
     * destructively parses a data string with the form "[int traceID, double impact]" to the traceTuple data
     * @param data the data string to be parsed
     */
    public void fromString(String data){
        String stringData = data.replace(" ", "");
        stringData = stringData.replace("[", "");
        stringData = stringData.replace("]", "");

        _traceID = Integer.parseInt(stringData.split(",")[0]);
        _impact = Double.parseDouble(stringData.split(",")[1]);
    }
}
