package proteinidentification;

import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import weka.classifiers.Classifier;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;

public class RALearningModel {

    Classifier classifier;
    Instances instances;

    public RALearningModel(Classifier c) {
        classifier = c;
    }

    public void buildModel(HashSet<RAProtein> positiveSet, HashSet<RAProtein> negativeSet) {
        try {
            FastVector attributes = new FastVector();
            FastVector classLabel = new FastVector();
            classLabel.addElement("Positive");
            classLabel.addElement("Negative");
            attributes.addElement(new Attribute("ClassLable", classLabel));
            int n = RAFeature.FEATURE_NUMBER-1;
            for (int i = 0; i < n; i++) {
                attributes.addElement(new Attribute("attributes_" + i));
            }
            instances = new Instances("ProteinFeatures", attributes, 0);
            instances.setClassIndex(0);
            addAttributes(classLabel.indexOf("Positive"), positiveSet);
            addAttributes(classLabel.indexOf("Negative"), negativeSet);
            classifier.buildClassifier(instances);
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }
    }

    public Hashtable<String, Double> applyModel(HashSet<RAProtein> proteinSet) {
        Hashtable<String, Double> result = new Hashtable<String, Double>();
        try {
            Iterator<RAProtein> iterator = proteinSet.iterator();
            while (iterator.hasNext()) {
                RAProtein protein = iterator.next();
                double[] attributes = new double[instances.numAttributes()];
                attributes[0] = Instance.missingValue();
                double[] features = protein.getFeature().getFeatureValuesBesidesCriterion();
                for (int i = 0; i < features.length; i++) {
                    attributes[i + 1] = features[i];
                }
                Instance instance = new Instance(1.0, attributes);
                instance.setDataset(instances);
                Double score = Double.valueOf(classifier.distributionForInstance(instance)[0]);
                result.put(protein.getProteinName(), score);
//                System.out.println(protein.getProteinName()+"\t"+score);
            }
        } catch (Exception e) {
            e.printStackTrace(System.err);
        } finally {
            return result;
        }
    }

    private void addAttributes(double classIndex, HashSet<RAProtein> proteinSet) {
        Iterator<RAProtein> iterator = proteinSet.iterator();
        while (iterator.hasNext()) {
            RAProtein protein = iterator.next();
//            System.out.println(classLabel+"\t"+protein.getProteinName());
            double[] attributes = new double[instances.numAttributes()];
            attributes[0] = classIndex;
            double[] features = protein.getFeature().getFeatureValuesBesidesCriterion();
            for (int i = 0; i < features.length; i++) {
                attributes[i + 1] = features[i];
            }
            instances.add(new Instance(1.0, attributes));
        }
    }
}
