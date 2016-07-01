//package proteinidentification;
//
//import java.util.ArrayList;
//import java.util.Collections;
//import java.util.Comparator;
//import java.util.Enumeration;
//import java.util.HashSet;
//import java.util.Hashtable;
//import java.util.List;
//import java.util.Map;
//
//public class RAEvaluator {
//
//    private HashSet<String> referenceSet;
//    private double[][] ROC;
//
//    public RAEvaluator(HashSet<String> rp) {
//        referenceSet = rp;
//    }
//
//    public void calculateROCCurvePoints(Hashtable<String, Double> data) {
////        System.out.println("************************************\n"+data);
//        HashSet<String> positiveSet = new HashSet<String>();
//        HashSet<String> negativeSet = new HashSet<String>();
//        Enumeration<String> e = data.keys();
//        while(e.hasMoreElements()){
//            String proteinName = e.nextElement();
//            if(referenceSet.contains(proteinName)){
//                positiveSet.add(proteinName);
//            }else{
//                negativeSet.add(proteinName);
//            }
//        }
//        ROC = new double[data.size()+1][2];
//        ROC[0][0] = 0.0;
//        ROC[0][1] = 0.0;
//        List<Map.Entry<String, Double>> list = null;
//        list = new ArrayList<Map.Entry<String, Double>>(data.entrySet());
//        HashtableComparator vc = new HashtableComparator();
//        Collections.sort(list, vc);
//
//        int i = 0, number = data.size(), n = 1;
//        int tp = 0,  fp = 0;
//        double tpr = 0, fpr = 0;
//        if (positiveSet.contains(list.get(0).getKey())) {
//            System.out.println(list.get(0).getValue());
//            tp++;
//        } else {
//            fp++;
//        }
//        tpr = tp * 1.0 / positiveSet.size();
//        fpr = fp * 1.0 / negativeSet.size();
//        ROC[1][0] = fpr;
//        ROC[1][1] = tpr;
//        double temp = list.get(0).getValue();
//        int cnt = 0;
//        for (i = 1; i < number; i++) {
//            if (positiveSet.contains(list.get(i).getKey())) {
//                if(list.get(i).getValue() == temp){
//                    cnt++;
//                }else{
//                    tp += (cnt + 1);
//                    cnt = 0;
//                }
//            } else {
//                fp++;
//            }
//            temp = list.get(i).getValue();
//            tpr = tp * 1.0 / positiveSet.size();
//            fpr = fp * 1.0 / negativeSet.size();
//            ROC[i + 1][0] = fpr;
//            ROC[i + 1][1] = tpr;
//        }
//    }
//    
//    public void exportROCCurvePoints(String output){
//        int i = 0, j = 0, k = 0, len = ROC.length;
//        double mx, my;
//        for (i = 1; i < len; i++) {
//            k = i;
//            mx = ROC[i][0];
//            my = ROC[i][1];
//            for (j = i; j > 0; j--) {
//                if (mx < ROC[j - 1][0]) {
//                    ROC[j][0] = ROC[j - 1][0];
//                    ROC[j][1] = ROC[j - 1][1];
//                    k = j - 1;
//                } else if (mx == ROC[j - 1][0]) {
//                    if (my < ROC[j - 1][1]) {
//                        ROC[j][0] = ROC[j - 1][0];
//                        ROC[j][1] = ROC[j - 1][1];
//                        k = j - 1;
//                    } else {
//                        break;
//                    }
//                } else {
//                    break;
//                }
//            }
//            if (k != i) {
//                ROC[k][0] = mx;
//                ROC[k][1] = my;
//            }
//        }
//        RAFileManager.writeROCToFile(output, ROC);
//    }
//    
//    class HashtableComparator implements Comparator<Map.Entry<String, Double>> {
//
//        @Override
//        public int compare(Map.Entry<String, Double> e1, Map.Entry<String, Double> e2) {
//            if (e2.getValue().doubleValue() > e1.getValue().doubleValue()) {
//                return 1;
//            } else if (e2.getValue().doubleValue() == e1.getValue().doubleValue()) {
//                return 0;
//            } else {
//                return -1;
//            }
//        }
//    }
//}
package proteinidentification;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;

public class RAEvaluator {

    private HashSet<String> referenceSet;
    private List<Roc> ROC;

    public RAEvaluator(HashSet<String> rp) {
        referenceSet = rp;
    }

    public void calculateROCCurvePoints(Hashtable<String, Double> data) {
//        System.out.println("************************************\n"+data);
        HashSet<String> positiveSet = new HashSet<String>();
        HashSet<String> negativeSet = new HashSet<String>();
        Enumeration<String> e = data.keys();
        while(e.hasMoreElements()){
            String proteinName = e.nextElement();
            if(referenceSet.contains(proteinName)){
                positiveSet.add(proteinName);
            }else{
                negativeSet.add(proteinName);
            }
        }
        
        List<Map.Entry<String, Double>> list = null;
        list = new ArrayList<Map.Entry<String, Double>>(data.entrySet());
        HashtableComparator vc = new HashtableComparator();
        Collections.sort(list, vc);

        int i = 0, number = data.size(), n = 1;
        int tp = 0,  fp = 0;
        double tpr = 0, fpr = 0;
        double temp = -1;
        ROC=new ArrayList<Roc>();
        for (i = 0; i < number; i++) {
            if(temp!=list.get(i).getValue()){
                if (positiveSet.contains(list.get(i).getKey())) {
                    tp++;
                } else {
                    fp++;
                }
                tpr = tp * 1.0 / positiveSet.size();
                fpr = fp * 1.0 / negativeSet.size();
                Roc roc=new Roc(fpr,tpr);
                ROC.add(roc);
            }
            else{
                    if (positiveSet.contains(list.get(i).getKey())) {
                        tp++;
                    } else {
                        fp++;
                    }
            }
            temp=list.get(i).getValue();
        }
        tpr = tp * 1.0 / positiveSet.size();
        fpr = fp * 1.0 / negativeSet.size();
        Roc roc=new Roc(fpr,tpr);
        ROC.add(roc);
    }
    
    public void exportROCCurvePoints(String output){
//        int i = 0, j = 0, k = 0, len = ROC.size();
//        double mx, my;
//        for (i = 1; i < len; i++) {
//            k = i;
//            mx = ROC[i][0];
//            my = ROC[i][1];
//            for (j = i; j > 0; j--) {
//                if (mx < ROC[j - 1][0]) {
//                    ROC[j][0] = ROC[j - 1][0];
//                    ROC[j][1] = ROC[j - 1][1];
//                    k = j - 1;
//                } else if (mx == ROC[j - 1][0]) {
//                    if (my < ROC[j - 1][1]) {
//                        ROC[j][0] = ROC[j - 1][0];
//                        ROC[j][1] = ROC[j - 1][1];
//                        k = j - 1;
//                    } else {
//                        break;
//                    }
//                } else {
//                    break;
//                }
//            }
//            if (k != i) {
//                ROC[k][0] = mx;
//                ROC[k][1] = my;
//            }
//        }
//        
//        
        RAFileManager.writeROCToFile(output, ROC);
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
