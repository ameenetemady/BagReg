
package proteinidentification;

public class RAProtein {
    private String proteinName;
    private RAFeature feature;
    
    public RAProtein(String name){
        proteinName=name;
        feature=RAFileManager.getFeature(name);
        
//        System.out.print("\t\""+name+"\"\t");
//        feature.printFeatureValues();
    }
    public String getProteinName(){
        return proteinName;
    }
    public RAFeature getFeature(){
        return feature;
    }
}
