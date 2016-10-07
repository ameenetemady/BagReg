package proteinidentification;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;
import weka.classifiers.Classifier;
import weka.classifiers.bayes.BayesNet;
import weka.classifiers.functions.Logistic;

public class RAAlgorithm {

    Classifier classifier;
    RALearningModel model;
    private HashSet<String> proteinSet;
    private HashSet<RAProtein> inferenceSet, positiveSet, negativeSet;
    private Hashtable<String, Double> result;
    private final double POSITIVE_RATIO = 0.2,NEGATIVE_RATIO = 0.2;

    public RAAlgorithm(Classifier c) {
        inferenceSet = new HashSet<RAProtein>();
        classifier = c;
        result = new Hashtable<String, Double>();
    }

    public void input(String input) {
        RAFileManager.inputPath = input;
//        inferenceSet = new HashSet<RAProtein>();
        proteinSet = RAFileManager.getProteinSet();
//        System.out.println("******************ProteinSet Number:"+proteinSet.size());          //////////
        Iterator<String> iterator = proteinSet.iterator();
        while (iterator.hasNext()) {
            String proteinName = iterator.next();
            RAProtein protein = new RAProtein(proteinName);
            inferenceSet.add(protein);
//            break;                    //////////////////////////////////////
        }
    }

    public void process() {
        ArrayList<String> featureNames = RAFeature.getFeatureNames();
        int number = featureNames.size();
        for (int i = 0; i < number; i++) {
            String name = featureNames.get(i);
            System.out.println("---------"+name+"---------");               ///////////
            Iterator<RAProtein> iterator = inferenceSet.iterator();
            while (iterator.hasNext()) {
                RAProtein protein = iterator.next();
                protein.getFeature().setCriterion(name);
            }
                buildTrainingSet();
                model = new RALearningModel(classifier);
                model.buildModel(positiveSet, negativeSet);
                Hashtable<String, Double> tempResult = model.applyModel(inferenceSet);
                addToFinalResult(tempResult);
//                RAFileManager.writeResultToFile("exportdata/"+name+".csv", new ArrayList<Map.Entry<String,Double>>(tempResult.entrySet()));
        }
        Enumeration<String> e = result.keys();
        while (e.hasMoreElements()) {
            String proteinName = e.nextElement();
            Double score = result.get(proteinName);
            score = Double.valueOf(score.doubleValue()/6.0);
            result.put(proteinName, score);
        }
    }

    public void output(String output) {
        ArrayList<Map.Entry<String, Double>> list = new ArrayList<Map.Entry<String, Double>>(result.entrySet());
        HashtableComparator comparator = new HashtableComparator();
        Collections.sort(list, comparator);
        RAFileManager.writeResultToFile(output, list);
    }

    private void buildTrainingSet() {
        positiveSet = new HashSet<RAProtein>();
        negativeSet = new HashSet<RAProtein>();
        ArrayList<RAProtein> list = new ArrayList<RAProtein>(inferenceSet);
        ProteinSetComparator comparator = new ProteinSetComparator();
        Collections.sort(list, comparator);
        int amount = list.size(), positiveNumber = (int) (POSITIVE_RATIO * amount),
                negativeNumber = (int)(NEGATIVE_RATIO*amount);
//        System.out.println("amount:"+amount+"\tnumber"+number);
        for (int i = 0; i < positiveNumber; i++) {
//            System.out.println("positive:"+list.get(i).getProteinName() +"("+list.get(i).getFeature().getCriterion()
//                    +")\tnegative:"+list.get(amount-1-i).getProteinName()+"("+list.get(amount-1-i).getFeature().getCriterion()+")");
            positiveSet.add(list.get(i));
        }
        for(int i=0;i<negativeNumber;i++){
            negativeSet.add(list.get(amount - 1 - i));
        }
        int p=positiveNumber,n=list.size()-negativeNumber-1;
        while(p<n && list.get(p-1).getFeature().getCriterion()==list.get(p).getFeature().getCriterion()){
            positiveSet.add(list.get(p));
            p++;
        }
        while(n>p && list.get(n+1).getFeature().getCriterion()==list.get(n).getFeature().getCriterion()){
            negativeSet.add(list.get(n));
            n--;
        }
//        printPositiveAndNegativeSet();
    }

    private void printPositiveAndNegativeSet(){
        Iterator<RAProtein> iterator = positiveSet.iterator();
        System.out.println("PositiveSet:");
        while(iterator.hasNext()){
            RAProtein protein = iterator.next();
            System.out.println("\t"+protein.getProteinName()+"\t("+protein.getFeature().getCriterion()+")");       //////
        }
        iterator = negativeSet.iterator();
        System.out.println("NegativeSet:");
        while(iterator.hasNext()){
            RAProtein protein = iterator.next();
            System.out.println("\t"+protein.getProteinName()+"\t("+protein.getFeature().getCriterion()+")");       //////
        }
        
    }
    private void addToFinalResult(Hashtable<String, Double> tempResult) {
        Enumeration<String> e = tempResult.keys();
        while (e.hasMoreElements()) {
            String proteinName = e.nextElement();
            Double finalScore = Double.valueOf(0);
            if (result.containsKey(proteinName)) {
                finalScore = result.get(proteinName);
            }
            finalScore = Double.valueOf(finalScore.doubleValue()+tempResult.get(proteinName).doubleValue());
            result.put(proteinName, finalScore);
        }
    }

    public static void main(String[] args) {
/*
        String[] paths = new String[6];
        paths[0] = "importdata/18_mixtures.txt" + " " + "exportdata/18_mixtures_PIR3A.csv" + " "
                + "importdata/18_mixtures_reference.csv" + " "
                + "importdata/18_mixtures_Prophet.csv" + " " + "exportdata/18_mixtures_Prophet_ROC_result.csv" + " "
                + "importdata/18_mixtures_MSBayesPro.csv" + " " + "exportdata/18_mixtures_MSBayesPro_ROC_result.csv" + " "
                + "importdata/18_mixtures_Fido.csv" + " " + "exportdata/18_mixtures_Fido_ROC_result.csv" + " "
                + "importdata/18_mixtures_PLP.csv" + " " + "exportdata/18_mixtures_PLP_ROC_result.csv" + " "
                + "exportdata/18_mixtures_PIR3A.csv" + " " + "exportdata/18_mixtures_PIR3A_ROC_result.csv";

        paths[1] = "importdata/DME.txt" + " " + "exportdata/DME_PIR3A.csv" + " "
                + "importdata/DME_reference.csv" + " "
                + "importdata/DME_Prophet.csv" + " " + "exportdata/DME_Prophet_ROC_result.csv" + " "
                + "importdata/DME_MSBayesPro.csv" + " " + "exportdata/DME_MSBayesPro_ROC_result.csv" + " "
                + "importdata/DME_Fido.csv" + " " + "exportdata/DME_Fido_ROC_result.csv" + " "
                + "importdata/DME_PLP.csv" + " " + "exportdata/DME_PLP_ROC_result.csv" + " "
                + "exportdata/DME_PIR3A.csv" + " " + "exportdata/DME_PIR3A_ROC_result.csv";

        paths[2] = "importdata/HumanEKC.txt" + " " + "exportdata/HumanEKC_PIR3A.csv" + " "
                + "importdata/HumanEKC_reference.csv" + " "
                + "importdata/HumanEKC_Prophet.csv" + " " + "exportdata/HumanEKC_Prophet_ROC_result.csv" + " "
                + "importdata/HumanEKC_MSBayesPro.csv" + " " + "exportdata/HumanEKC_MSBayesPro_ROC_result.csv" + " "
                + "importdata/HumanEKC_Fido.csv" + " " + "exportdata/HumanEKC_Fido_ROC_result.csv" + " "
                + "importdata/HumanEKC_PLP.csv" + " " + "exportdata/HumanEKC_PLP_ROC_result.csv" + " "
                + "exportdata/HumanEKC_PIR3A.csv" + " " + "exportdata/HumanEKC_PIR3A_ROC_result.csv";

        paths[3] = "importdata/HumanMD.txt" + " " + "exportdata/HumanMD_PIR3A.csv" + " "
                + "importdata/HumanMD_reference.csv" + " "
                + "importdata/HumanMD_Prophet.csv" + " " + "exportdata/HumanMD_Prophet_ROC_result.csv" + " "
                + "importdata/HumanMD_MSBayesPro.csv" + " " + "exportdata/HumanMD_MSBayesPro_ROC_result.csv" + " "
                + "importdata/HumanMD_Fido.csv" + " " + "exportdata/HumanMD_Fido_ROC_result.csv" + " "
                + "importdata/HumanMD_PLP.csv" + " " + "exportdata/HumanMD_PLP_ROC_result.csv" + " "
                + "exportdata/HumanMD_PIR3A.csv" + " " + "exportdata/HumanMD_PIR3A_ROC_result.csv";

        paths[4] = "importdata/Sigma_49.txt" + " " + "exportdata/Sigma_49_PIR3A.csv" + " "
                + "importdata/Sigma_49_reference.csv" + " "
                + "importdata/Sigma_49_Prophet.csv" + " " + "exportdata/Sigma_49_Prophet_ROC_result.csv" + " "
                + "importdata/Sigma_49_MSBayesPro.csv" + " " + "exportdata/Sigma_49_MSBayesPro_ROC_result.csv" + " "
                + "importdata/Sigma_49_Fido.csv" + " " + "exportdata/Sigma_49_Fido_ROC_result.csv" + " "
                + "importdata/Sigma_49_PLP.csv" + " " + "exportdata/Sigma_49_PLP_ROC_result.csv" + " "
                + "exportdata/Sigma_49_PIR3A.csv" + " " + "exportdata/Sigma_49_PIR3A_ROC_result.csv";

        paths[5] = "importdata/Yeast.txt" + " " + "exportdata/Yeast_PIR3A.csv" + " "
                + "importdata/Yeast_reference.csv" + " "
                + "importdata/Yeast_Prophet.csv" + " " + "exportdata/Yeast_Prophet_ROC_result.csv" + " "
                + "importdata/Yeast_MSBayesPro.csv" + " " + "exportdata/Yeast_MSBayesPro_ROC_result.csv" + " "
                + "importdata/Yeast_Fido.csv" + " " + "exportdata/Yeast_Fido_ROC_result.csv" + " "
                + "importdata/Yeast_PLP.csv" + " " + "exportdata/Yeast_PLP_ROC_result.csv" + " "
                + "exportdata/Yeast_PIR3A.csv" + " " + "exportdata/Yeast_PIR3A_ROC_result.csv";

        for (int i = 0; i < 6; i++) {
            String[] cut = paths[i].split(" ");
            System.out.println(cut[0]);
            RAAlgorithm algorithm = new RAAlgorithm(new BayesNet());    //change classifier algorithm
            algorithm.input(cut[0]);
            algorithm.process();
            algorithm.output(cut[1]);
            RAEvaluator evaluator = new RAEvaluator(RAFileManager.getReferenceSet(cut[2]));
            for (int j = 0; j < 5; j++) {
                evaluator.calculateROCCurvePoints(RAFileManager.getExistingResult(cut[3 + 2 * j]));
                evaluator.exportROCCurvePoints(cut[4 + 2 * j]);
            }
        }
	*/

	 String peptideFile=args[0];
	 String resultFile=args[1];

         RAAlgorithm algorithm = new RAAlgorithm(new BayesNet());    //change classifier algorithm
         algorithm.input(peptideFile);
         algorithm.process();
         algorithm.output(resultFile);
    }

    class ProteinSetComparator implements Comparator<RAProtein> {

        @Override
        public int compare(RAProtein p1, RAProtein p2) {
            if (p2.getFeature().getCriterion() > p1.getFeature().getCriterion()) {
                return 1;
            } else if (p2.getFeature().getCriterion() == p1.getFeature().getCriterion()) {
                return 0;
            } else {
                return -1;
            }
        }
    }

    class HashtableComparator implements Comparator<Map.Entry<String, Double>> {

        @Override
        public int compare(Map.Entry<String, Double> e1, Map.Entry<String, Double> e2) {
            if (e2.getValue().doubleValue() > e1.getValue().doubleValue()) {
                return 1;
            } else if (e2.getValue().doubleValue() == e1.getValue().doubleValue()) {
                return 0;
            } else {
                return -1;
            }
        }
    }
}
